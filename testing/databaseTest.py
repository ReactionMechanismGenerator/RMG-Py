#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2019 Prof. William H. Green (whgreen@mit.edu),           #
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

"""
This scripts runs tests on the database
"""

import itertools
import logging
from copy import copy
from collections import defaultdict

import nose
import nose.tools
import numpy as np
import quantities as pq

import rmgpy.kinetics
import rmgpy.constants
from rmgpy import settings
from rmgpy.data.base import LogicOr
from rmgpy.data.rmg import RMGDatabase
from rmgpy.exceptions import ImplicitBenzeneError, UnexpectedChargeError
from rmgpy.molecule import Group
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.molecule.pathfinder import find_shortest_path
from rmgpy.quantity import ScalarQuantity


class TestDatabase(object):  # cannot inherit from unittest.TestCase if we want to use nose test generators
    """
    Contains unit tests for the database for rigorous error checking.
    """

    @classmethod
    def setUpClass(cls):
        """
        Load the database before running the tests.
        """
        database_directory = settings['database.directory']
        cls.database = RMGDatabase()
        cls.database.load(database_directory, kinetics_families='all')

    # These are generators, that call the methods below.
    def test_kinetics(self):
        for family_name, family in self.database.kinetics.families.items():

            test = lambda x: self.kinetics_check_correct_number_of_nodes_in_rules(family_name)
            test_name = "Kinetics family {0}: rules have correct number of nodes?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, None

            test = lambda x: self.kinetics_check_nodes_in_rules_found_in_groups(family_name)
            test_name = "Kinetics family {0}: rules' nodes exist in the groups?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, None

            test = lambda x: self.kinetics_check_groups_found_in_tree(family_name)
            test_name = "Kinetics family {0}: groups are in the tree with proper parents?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, None

            test = lambda x: self.kinetics_check_groups_nonidentical(family_name)
            test_name = "Kinetics family {0}: groups are not identical?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_check_child_parent_relationships(family_name)
            test_name = "Kinetics family {0}: parent-child relationships are correct?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_check_siblings_for_parents(family_name)
            test_name = "Kinetics family {0}: sibling relationships are correct?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_check_cd_atom_type(family_name)
            test_name = "Kinetics family {0}: Cd, CS, CO, and Cdd atomtype used correctly?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_check_reactant_and_product_template(family_name)
            test_name = "Kinetics family {0}: reactant and product templates correctly defined?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            # tests for surface families
            if 'surface' in family_name.lower():
                test = lambda x: self.kinetics_check_surface_training_reactions_can_be_used(family_name)
                test_name = "Kinetics surface family {0}: entries can be used to generate rate rules?".format(family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

                test = lambda x: self.kinetics_check_training_reactions_have_surface_attributes(family_name)
                test_name = "Kinetics surface family {0}: entries have surface attributes?".format(family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

            # these families have some sort of difficulty which prevents us from testing accessibility right now
            difficult_families = ['Diels_alder_addition', 'Intra_R_Add_Exocyclic', 'Intra_R_Add_Endocyclic', 'Retroene']
            generated_trees = ["R_Recombination"]

            if len(family.forward_template.reactants) < len(family.groups.top) and family_name not in difficult_families:
                test = lambda x: self.kinetics_check_unimolecular_groups(family_name)
                test_name = "Kinetics family {0} check that unimolecular group is formatted correctly?".format(
                    family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

            if family_name not in difficult_families and family_name not in generated_trees:
                test = lambda x: self.kinetics_check_sample_descends_to_group(family_name)
                test_name = "Kinetics family {0}: Entry is accessible?".format(family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

                test = lambda x: self.kinetics_check_sample_can_react(family_name)
                test_name = "Kinetics family {0}: Recipe applies to group entry?".format(family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

            for depository in family.depositories:
                test = lambda x: self.kinetics_check_adjlists_nonidentical(depository)
                test_name = "Kinetics depository {0}: check adjacency lists are nonidentical?".format(depository.label)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, depository.label

                test = lambda x: self.kinetics_check_rate_units_are_correct(depository, tag='depository')
                test_name = "Kinetics depository {0}: check rates have correct units?".format(depository.label)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, depository.label

        for library_name, library in self.database.kinetics.libraries.items():
            test = lambda x: self.kinetics_check_adjlists_nonidentical(library)
            test_name = "Kinetics library {0}: check adjacency lists are nonidentical?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

            test = lambda x: self.kinetics_check_rate_units_are_correct(library)
            test_name = "Kinetics library {0}: check rates have correct units?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

            test = lambda x: self.kinetics_check_library_rates_are_reasonable(library)
            test_name = "Kinetics library {0}: check rates are reasonable?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

            # tests for surface families
            if 'surface' in library_name.lower():
                test = lambda x: self.kinetics_check_surface_library_reactions_have_surface_attributes(library)
                test_name = "Kinetics surface library {0}: entries have surface attributes?".format(library_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

    def test_thermo(self):
        for group_name, group in self.database.thermo.groups.items():
            test = lambda x: self.general_check_nodes_found_in_tree(group_name, group)
            test_name = "Thermo groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_groups_nonidentical(group_name, group)
            test_name = "Thermo groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_child_parent_relationships(group_name, group)
            test_name = "Thermo groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_siblings_for_parents(group_name, group)
            test_name = "Thermo groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_cd_atom_type(group_name, group)
            test_name = "Thermo groups {0}: Cd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_sample_descends_to_group(group_name, group)
            test_name = "Thermo groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            # tests for adsorption groups
            if 'adsorption' in group_name.lower():
                test = lambda x: self.check_surface_thermo_groups_have_surface_attributes(group_name, group)
                test_name = "Thermo surface groups {0}: Entry has metal attributes?".format(group_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, group_name

        for library_name, library in self.database.thermo.libraries.items():
            if 'surface' in library_name.lower():
                test = lambda x: self.check_surface_thermo_libraries_have_surface_attributes(library_name, library)
                test_name = "Thermo surface libraries {0}: Entry has metal attributes?".format(library_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, group_name

    def test_solvation(self):
        for group_name, group in self.database.solvation.groups.items():
            test = lambda x: self.general_check_nodes_found_in_tree(group_name, group)
            test_name = "Solvation groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_groups_nonidentical(group_name, group)
            test_name = "Solvation groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_child_parent_relationships(group_name, group)
            test_name = "Solvation groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_siblings_for_parents(group_name, group)
            test_name = "Solvation groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_cd_atom_type(group_name, group)
            test_name = "Solvation groups {0}: Cd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_sample_descends_to_group(group_name, group)
            test_name = "Solvation groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

    def test_statmech(self):
        for group_name, group in self.database.statmech.groups.items():
            test = lambda x: self.general_check_nodes_found_in_tree(group_name, group)
            test_name = "Statmech groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_groups_nonidentical(group_name, group)
            test_name = "Statmech groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_child_parent_relationships(group_name, group)
            test_name = "Statmech groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_siblings_for_parents(group_name, group)
            test_name = "Statmech groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_cd_atom_type(group_name, group)
            test_name = "Statmech groups {0}: Cd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_sample_descends_to_group(group_name, group)
            test_name = "Statmech groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

    def test_transport(self):
        for group_name, group in self.database.transport.groups.items():
            test = lambda x: self.general_check_nodes_found_in_tree(group_name, group)
            test_name = "Transport groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_groups_nonidentical(group_name, group)
            test_name = "Transport groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_child_parent_relationships(group_name, group)
            test_name = "Transport groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_siblings_for_parents(group_name, group)
            test_name = "Transport groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_cd_atom_type(group_name, group)
            test_name = "Transport groups {0}: Cd, CS, CO, and Cdd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_check_sample_descends_to_group(group_name, group)
            test_name = "Transport groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

    def test_metal_libraries(self):
        for library_name, library in self.database.thermo.surface['metal'].libraries.items():
            test = lambda x: self.general_check_metal_database_has_catalyst_properties(library)
            test_name = "Metal library {0}: Entries have catalyst properties?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

            test = lambda x: self.general_check_metal_database_has_reasonable_labels(library)
            test_name = "Metal library {0}: Entries have reasonable labels?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

    # These are the actual tests, that don't start with a "test_" name:
    def kinetics_check_surface_training_reactions_can_be_used(self, family_name):
        """Test that surface training reactions can be averaged and used for generating rate rules"""
        family = self.database.kinetics.families[family_name]
        family.add_rules_from_training(thermo_database=self.database.thermo)

    def general_check_metal_database_has_catalyst_properties(self, library):
        """Test that each entry has catalyst properties"""
        for entry in library.entries.values():
            if not entry.binding_energies:
                raise AttributeError('Entry {} has no binding energies'.format(entry.label))
            assert isinstance(entry.binding_energies, dict)
            for element in 'CHON':
                if not entry.binding_energies[element]:
                    raise KeyError('Entry {} has no {} binding energy'.format(entry.label, element))
                if not isinstance(entry.binding_energies[element], ScalarQuantity):
                    raise TypeError('Entry {} binding energy value for {} should be a ScalarQuantity, but is type {}'.format(
                        entry.label, element, type(entry.binding_energies[element])))
                if not isinstance(entry.binding_energies[element].value, float):
                    raise TypeError('Entry {} binding energy for {} should be a float, but is type {}'.format(
                        entry.label, element, type(entry.binding_energies[element].value)))
                assert entry.binding_energies[element].value < 0.  # binding energies should all be negative... probably
                assert entry.binding_energies[element].units == 'eV/molecule'

            if not entry.surface_site_density:
                raise AttributeError('Entry {} has no surface site density'.format(entry.label))
            assert isinstance(entry.surface_site_density, ScalarQuantity)
            if not isinstance(entry.surface_site_density.value, float):
                raise TypeError('Entry {} should be a float, but is type {}'.format(entry.label, type(entry.surface_site_density.value)))
            if not isinstance(entry.surface_site_density.units, str):
                raise TypeError('Entry {} should be a str, but is type {}'.format(entry.label, type(entry.surface_site_density.units)))
            assert 1e-4 > entry.surface_site_density.value_si > 1e-6  # values should be reasonable

            assert isinstance(entry.metal, str)  # all entries should have a metal attribute, at minimum
            if entry.facet:
                assert isinstance(entry.facet, str)
            if entry.site:
                assert isinstance(entry.site, str)

    def general_check_metal_database_has_reasonable_labels(self, library):
        """Test that each entry has a reasonable label corresponding to its metal and facet"""
        for entry in library.entries.values():
            if entry.metal not in entry.label:
                raise NameError('Entry {} with metal attribute {} does not have metal in its label'.format(entry.label, entry.metal))
            if entry.facet not in entry.label:
                raise NameError('Entry {} with facet attribute {} does not have facet in its label'.format(entry.label, entry.facet))
            if not entry.label[0].isupper():
                raise NameError('Entry {} should start with a capital letter'.format(entry.label))

    def kinetics_check_training_reactions_have_surface_attributes(self, family_name):
        """Test that each surface training reaction has surface attributes"""
        family = self.database.kinetics.families[family_name]
        training = family.get_training_depository().entries.values()
        failed = False
        for entry in training:
            if not entry.metal:
                logging.error(f'Expected a metal attribute for {entry} in {family} family but found {entry.metal!r}')
                failed = True
            else:
                assert isinstance(entry.metal, str)

            if entry.facet:
                assert isinstance(entry.facet, str)
            if entry.site:
                assert isinstance(entry.site, str)

        if failed:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_surface_library_reactions_have_surface_attributes(self, library):
        """Test that each surface reaction library has surface attributes"""
        entries = library.entries.values()
        failed = False
        if '_Pt' in library.label:
            for entry in entries:
                if entry.metal is not 'Pt':
                    logging.error(f'Expected {entry} metal attribute in {library} library to match Pt, but was {entry.metal}')
                    failed = True
        if '_Ni' in library.label:
            for entry in entries:
                if entry.metal is not 'Ni':
                    logging.error(f'Expected {entry} metal attribute in {library} library to match Ni, but was {entry.metal}')
                    failed = True
        for entry in entries:
            if isinstance(entry.metal, type(None)):
                logging.error(f'Expected a metal attribute in {library} library for {entry} but found None')
                failed = True
        if failed:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_correct_number_of_nodes_in_rules(self, family_name):
        """
        This test ensures that each rate rule contains the proper number of nodes according to the family it originates.
        """
        family = self.database.kinetics.families[family_name]
        expected_number_nodes = len(family.get_root_template())
        tst = []
        for label, entries in family.rules.entries.items():
            for entry in entries:
                nodes = label.split(';')
                tst.append((len(nodes), expected_number_nodes,
                            "Wrong number of groups or semicolons in family {family} rule {entry}. Should be "
                            "{num_nodes}".format(family=family_name, entry=entry, num_nodes=expected_number_nodes)))

        boo = False
        for item in tst:
            if item[0] != item[1]:
                boo = True
                logging.error(item[2])
        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_nodes_in_rules_found_in_groups(self, family_name):
        """
        This test ensures that each rate rule contains nodes that exist in the
        groups and that they match the order of the forwardTemplate.
        """
        family = self.database.kinetics.families[family_name]

        # List of the each top node's descendants (including the top node)
        top_descendants = []
        for topNode in family.get_root_template():
            nodes = [topNode]
            nodes.extend(family.groups.descendants(topNode))
            top_descendants.append(nodes)

        top_group_order = ';'.join(topNode.label for topNode in family.get_root_template())
        tst1 = []
        tst2 = []
        for label, entries in family.rules.entries.items():
            for entry in entries:
                nodes = label.split(';')
                for i, node in enumerate(nodes):
                    tst1.append((node in family.groups.entries,
                                 "In {family} family, no group definition found for label {label} in rule "
                                 "{entry}".format(family=family_name, label=node, entry=entry)))
                    tst2.append((family.groups.entries[node] in top_descendants[i],
                                 "In {family} family, rule {entry} was found with groups out of order. "
                                 "The correct order for a rule should be subgroups of {top}.".format(
                                     family=family_name, entry=entry, top=top_group_order)))
        boo = False
        for i in range(len(tst1)):
            if not tst1[i][0]:
                logging.error(tst1[i][1])
                boo = True
            if not tst2[i][0]:
                logging.error(tst2[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_groups_found_in_tree(self, family_name):
        """
        This test checks whether groups are found in the tree, with proper parents.
        """
        family = self.database.kinetics.families[family_name]
        tst = []
        tst1 = []
        tst2 = []
        tst3 = []
        for nodeName, nodeGroup in family.groups.entries.items():
            tst.append(('[' in nodeName or ']' in nodeName,
                        "Group {group} in {family} family contains square brackets [ ] in the label, which are "
                        "not allowed.".format(group=nodeName, family=family_name)))
            ascend_parent = nodeGroup

            # Check whether the node has proper parents unless it is the top reactant or product node
            while ascend_parent not in family.groups.top and ascend_parent not in family.forward_template.products:
                child = ascend_parent
                ascend_parent = ascend_parent.parent
                tst1.append((ascend_parent is not None,
                             "Group {group} in {family} family was found in the tree without a proper parent.".format(
                                 group=child, family=family_name)))
                tst2.append((child in ascend_parent.children,
                             "Group {group} in {family} family was found in the tree without a proper parent.".format(
                                 group=nodeName, family=family_name)))
                tst3.append((child is ascend_parent,
                             "Group {group} in {family} family is a parent to itself".format(
                                 group=nodeName, family=family_name)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True
        for i in range(len(tst1)):
            if not tst1[i][0]:
                logging.error(tst1[i][1])
                boo = True
            if not tst2[i][0]:
                logging.error(tst2[i][1])
                boo = True
            if tst3[i][0]:
                logging.error(tst3[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_groups_nonidentical(self, family_name):
        """
        This test checks that the groups are non-identical.
        """
        from rmgpy.data.base import Database
        original_family = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = original_family.groups.entries
        entries_copy = copy(family.entries)
        tst = []
        for nodeName, nodeGroup in family.entries.items():
            del entries_copy[nodeName]
            for nodeNameOther, nodeGroupOther in entries_copy.items():
                tst.append((family.match_node_to_node(nodeGroup, nodeGroupOther),
                            "Group {group} in {family} family was found to be identical to group {groupOther}".format(
                                group=nodeName, family=family_name, groupOther=nodeNameOther)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_child_parent_relationships(self, family_name):
        """
        This test checks that groups' parent-child relationships are correct in the database.
        """
        from rmgpy.data.base import Database
        original_family = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = original_family.groups.entries
        tst = []
        for nodeName, childNode in family.entries.items():
            # top nodes and product nodes don't have parents by definition, so they get an automatic pass:
            if childNode in original_family.groups.top or childNode in original_family.forward_template.products:
                continue
            parent_node = childNode.parent

            if parent_node is None:
                # This is a mistake in the database, but it should be caught by kinetics_checkGroupsFoundInTree
                # so rather than report it twice or crash, we'll just silently carry on to the next node.
                continue
            elif parent_node in original_family.forward_template.products:
                # This is a product node made by training reactions which we do not need to heck
                continue
            # Check whether the node has proper parents unless it is the top reactant or product node
            # The parent should be more general than the child
            tst.append((family.match_node_to_child(parent_node, childNode),
                        "In {family} family, group {parent} is not a proper parent of its child "
                        "{child}.".format(family=family_name, parent=parent_node, child=nodeName)))
            # check that parentNodes which are LogicOr do not have an ancestor that is a Group
            # If it does, then the child_node must also be a child of the ancestor
            if isinstance(parent_node.item, LogicOr):
                ancestor_node = parent_node
                while ancestor_node not in original_family.groups.top and isinstance(ancestor_node.item, LogicOr):
                    ancestor_node = ancestor_node.parent
                if isinstance(ancestor_node.item, Group):
                    tst.append((family.match_node_to_child(ancestor_node, childNode),
                                "In {family} family, group {ancestor} is not a proper ancestor of its child "
                                "{child}.".format(family=family_name, ancestor=ancestor_node, child=nodeName)))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_siblings_for_parents(self, family_name):
        """
        This test checks that siblings in a tree are not actually parent/child

        See general_checkSiblingsForParents comments for more detailed description
        of the test.
        """
        from rmgpy.data.base import Database
        original_family = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = original_family.groups.entries
        tst = []
        for nodeName, node in family.entries.items():
            # Some families also construct a 2-level trees for the products
            # (root with all entries down one level) We don't care about this
            # tree as it is not used in searching, so we ignore products
            if node in original_family.forward_template.products:
                continue
            for index, child1 in enumerate(node.children):
                for child2 in node.children[index + 1:]:
                    tst.append((family.match_node_to_child(child1, child2),
                                "In family {0}, node {1} is a parent of {2}, but they are written as siblings.".format(
                                    family_name, child1, child2)))
        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_adjlists_nonidentical(self, database):
        """
        This test checks whether adjacency lists of reactants in a KineticsDepository or KineticsLibrary
        database object are nonidentical.
        """
        species_dict = {}
        entries = database.entries.values()
        for entry in entries:
            for reactant in entry.item.reactants:
                if reactant.label not in species_dict:
                    species_dict[reactant.label] = reactant

            for product in entry.item.products:
                if product.label not in species_dict:
                    species_dict[product.label] = product

        tst = []
        boo = False
        # Go through all species to make sure they are nonidentical
        species_list = list(species_dict.values())
        labeled_atoms = [species.molecule[0].get_all_labeled_atoms() for species in species_list]
        for i in range(len(species_list)):
            for j in range(i + 1, len(species_list)):
                initial_map = {}
                try:
                    atom_labels = set(list(labeled_atoms[i].keys()) +
                                      list(labeled_atoms[j].keys()))
                    for atomLabel in atom_labels:
                        initial_map[labeled_atoms[i][atomLabel]] = labeled_atoms[j][atomLabel]
                except KeyError:
                    # atom labels did not match, therefore not a match
                    continue
                m1 = species_list[i].molecule[0]
                m2 = species_list[j].molecule[0]
                if not m1.is_mapping_valid(m2, initial_map, equivalent=True):
                    # the mapping is invalid so they're not isomorphic
                    continue
                if m1.is_isomorphic(m2, initial_map):
                    logging.error("Species {0} and species {1} in {2} database were found to be identical.".format(
                                species_list[i].label, species_list[j].label, database.label))
                    boo = True
        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_rate_units_are_correct(self, database, tag='library'):
        """
        This test ensures that every reaction has acceptable units on the A factor.
        """
        boo = False

        dimensionalities = {
            1: (1 / pq.s).dimensionality,
            2: (pq.m ** 3 / pq.mole / pq.s).dimensionality,
            3: ((pq.m ** 6) / (pq.mole ** 2) / pq.s).dimensionality,
        }

        for entry in database.entries.values():
            k = entry.data
            rxn = entry.item
            molecularity = len(rxn.reactants)
            surface_reactants = sum([1 for s in rxn.reactants if s.contains_surface_site()])
            try:

                if isinstance(k, rmgpy.kinetics.StickingCoefficient):
                    "Should be dimensionless"
                    a_factor = k.A
                    if a_factor.units:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid units {3}'.format(
                            rxn, tag, database.label, a_factor.units))
                elif isinstance(k, rmgpy.kinetics.SurfaceArrhenius):
                    a_factor = k.A
                    expected = copy(dimensionalities[molecularity])
                    # for each surface reactant but one, switch from (m3/mol) to (m2/mol)
                    for _ in range(surface_reactants - 1):
                        expected[pq.m] -= 1
                    if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != expected:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid units {3}'.format(
                            rxn, tag, database.label, a_factor.units))
                elif isinstance(k, rmgpy.kinetics.Arrhenius):  # (but not SurfaceArrhenius, which came first)
                    a_factor = k.A
                    if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != dimensionalities[molecularity]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid units {3}'.format(
                            rxn, tag, database.label, a_factor.units))
                elif isinstance(k, (rmgpy.kinetics.Lindemann, rmgpy.kinetics.Troe)):
                    a_factor = k.arrheniusHigh.A
                    if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != dimensionalities[molecularity]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid high-pressure limit units {3}'.format(
                            rxn, tag, database.label, a_factor.units))
                elif isinstance(k, (rmgpy.kinetics.Lindemann, rmgpy.kinetics.Troe, rmgpy.kinetics.ThirdBody)):
                    a_factor = k.arrheniusLow.A
                    if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != dimensionalities[molecularity + 1]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid low-pressure limit units {3}'.format(
                            rxn, tag, database.label, a_factor.units))
                elif hasattr(k, 'highPlimit') and k.highPlimit is not None:
                    a_factor = k.highPlimit.A
                    if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != dimensionalities[molecularity - 1]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid high-pressure limit units {3}'.format(
                            rxn, tag, database.label, a_factor.units))
                elif isinstance(k, rmgpy.kinetics.MultiArrhenius):
                    for num, arrhenius in enumerate(k.arrhenius):
                        a_factor = arrhenius.A
                        if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != dimensionalities[molecularity]:
                            boo = True
                            logging.error(
                                'Reaction {0} from {1} {2}, has invalid units {3} on rate expression {4}'.format(
                                    rxn, tag, database.label, a_factor.units, num + 1)
                            )

                elif isinstance(k, rmgpy.kinetics.PDepArrhenius):
                    for pa, arrhenius in zip(k.pressures.value_si, k.arrhenius):
                        P = rmgpy.quantity.Pressure(1, k.pressures.units)
                        P.value_si = pa

                        if isinstance(arrhenius, rmgpy.kinetics.MultiArrhenius):
                            # A PDepArrhenius may have MultiArrhenius within it
                            # which is distinct (somehow) from MultiPDepArrhenius
                            for num, arrhenius2 in enumerate(arrhenius.arrhenius):
                                a_factor = arrhenius2.A
                                if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != \
                                        dimensionalities[molecularity]:
                                    boo = True
                                    logging.error(
                                        'Reaction {0} from {1} {2}, has invalid units {3} on {4!r} rate expression '
                                        '{5}'.format(rxn, tag, database.label, a_factor.units, P, num + 1)
                                    )
                        else:
                            a_factor = arrhenius.A
                            if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != \
                                    dimensionalities[molecularity]:
                                boo = True
                                logging.error(
                                    'Reaction {0} from {1} {2}, has invalid {3!r} units {4}'.format(
                                        rxn, tag, database.label, P, a_factor.units)
                                )

                elif isinstance(k, rmgpy.kinetics.MultiPDepArrhenius):
                    for num, k2 in enumerate(k.arrhenius):
                        for pa, arrhenius in zip(k2.pressures.value_si, k2.arrhenius):
                            P = rmgpy.quantity.Pressure(1, k2.pressures.units)
                            P.value_si = pa
                            if isinstance(arrhenius, rmgpy.kinetics.MultiArrhenius):
                                # A MultiPDepArrhenius may have MultiArrhenius within it
                                for arrhenius2 in arrhenius.arrhenius:
                                    a_factor = arrhenius2.A
                                    if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != \
                                            dimensionalities[molecularity]:
                                        boo = True
                                        logging.error(
                                            'Reaction {0} from {1} {2}, has invalid units {3} on {4!r} rate expression '
                                            '{5!r}'.format(rxn, tag, database.label, a_factor.units, P, arrhenius2)
                                        )
                            else:
                                a_factor = arrhenius.A
                                if pq.Quantity(1.0, a_factor.units).simplified.dimensionality != \
                                        dimensionalities[molecularity]:
                                    boo = True
                                    logging.error(
                                        'Reaction {0} from {1} {2}, has invalid {3!r} units {4} in rate expression '
                                        '{5}'.format(rxn, tag, database.label, P, a_factor.units, num)
                                    )

                elif isinstance(k, rmgpy.kinetics.Chebyshev):
                    if pq.Quantity(1.0, k.kunits).simplified.dimensionality != dimensionalities[molecularity]:
                        boo = True
                        logging.error(
                            'Reaction {0} from {1} {2}, has invalid units {3}'.format(
                                rxn, tag, database.label, k.kunits)
                        )

                else:
                    logging.warning('Reaction {0} from {1} {2}, did not have units checked.'.format(
                        rxn, tag, database.label))
            except:
                logging.error("Error when checking units on reaction {0} from {1} {2} with "
                              "rate expression {3!r}.".format(rxn, tag, database.label, k))
                raise
        if boo:
            raise ValueError('{0} {1} has some incorrect units'.format(tag.capitalize(), database.label))

    def kinetics_check_library_rates_are_reasonable(self, library):
        """
        This test ensures that every library reaction has reasonable kinetics at 1000 K, 1 bar
        """
        T = 1000.0  # K
        P = 100000.0  # 1 bar in Pa
        Na = rmgpy.constants.Na  # molecules/mol
        h_rad_mass = 1.6737236e-27  # kg
        h_rad_diam = 2.4e-10  # m
        kB = rmgpy.constants.kB  # m2 * kg *s^-2 * K^-1
        h = rmgpy.constants.h  # m2 kg / s
        boo = False
        tst_limit = (kB * T) / h
        collision_limit = Na * np.pi * h_rad_diam ** 2 * np.sqrt(8 * kB * T / (np.pi * h_rad_mass / 2))
        for entry in library.entries.values():
            if entry.item.is_surface_reaction():
                # Don't check surface reactions
                continue
            k = entry.data.get_rate_coefficient(T, P)
            rxn = entry.item
            if k < 0:
                boo = True
                logging.error('library reaction {0} from library {1}, had a negative rate at 1000 K, 1 bar'.format(
                    rxn, library.label))
            if len(rxn.reactants) == 1 and not rxn.allow_max_rate_violation:
                if k > tst_limit:
                    boo = True
                    logging.error('library reaction {0} from library {1}, exceeds the TST limit at 1000 K, 1 bar of '
                                  '{2} mol/(m3*s) at {3} mol/(m3*s) and did not have allow_max_rate_violation=True'
                                  ''.format(rxn, library.label, tst_limit, k))
            elif len(rxn.reactants) == 2 and not rxn.allow_max_rate_violation:
                if k > collision_limit:
                    boo = True
                    logging.error('library reaction {0} from library {1}, exceeds the collision limit at 1000 K, 1 bar '
                                  'of {2} mol/(m3*s) at {3} mol/(m3*s) and did not have allow_max_rate_violation=True'
                                  ''.format(rxn, library.label, collision_limit, k))
        if boo:
            raise ValueError('library {0} has unreasonable rates'.format(library.label))

    def kinetics_check_reactant_and_product_template(self, family_name):
        """
        This test checks whether the reactant and product templates within a family are correctly defined.
        For a reversible family, the reactant and product templates must have matching labels.
        For a non-reversible family, the reactant and product templates must have non-matching labels,
        otherwise overwriting may occur.
        """
        family = self.database.kinetics.families[family_name]
        if family.own_reverse:
            nose.tools.assert_equal(family.forward_template.reactants, family.forward_template.products)
        else:
            tst = []
            reactant_labels = [reactant.label for reactant in family.forward_template.reactants]
            product_labels = [product.label for product in family.forward_template.products]
            for reactant_label in reactant_labels:
                for product_label in product_labels:
                    tst.append((reactant_label == product_label,
                                "Reactant label {0} matches that of product label {1} in a non-reversible family "
                                "template.  Please rename product label.".format(reactant_label, product_label)))

            boo = False
            for i in range(len(tst)):
                if tst[i][0]:
                    logging.error(tst[i][1])
                    boo = True

            if boo:
                raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_cd_atom_type(self, family_name):
        """
        This test checks that groups containing Cd, CO, CS and Cdd atomtypes are used
        correctly according to their strict definitions
        """
        family = self.database.kinetics.families[family_name]
        target_label = ['Cd', 'CO', 'CS', 'Cdd']
        target_atom_types = [ATOMTYPES[x] for x in target_label]

        # ignore product entries that get created from training reactions
        ignore = []
        if not family.own_reverse:
            for product in family.forward_template.products:
                ignore.append(product)
                ignore.extend(product.children)
        else:
            ignore = []
        tst = []
        for entryName, entry in family.groups.entries.items():
            # ignore products
            if entry in ignore: continue
            # ignore LogicOr groups
            if isinstance(entry.item, Group):
                for index, atom in enumerate(entry.item.atoms):
                    for atomtype1 in atom.atomtype:
                        if atomtype1 in target_atom_types:
                            break
                    # If Cd not found in atomTypes, go to next atom
                    else:
                        continue
                    # Create list of all the atomTypes that should be present in addition or instead of Cd
                    correct_atom_list = []
                    num_of_d_bonds = sum(
                        [1 if x.order[0] is 'D' and len(x.order) == 1 else 0 for x in atom.bonds.values()])
                    if num_of_d_bonds == 2:
                        correct_atom_list.append('Cdd')
                    elif num_of_d_bonds == 1:
                        for ligand, bond in atom.bonds.items():
                            # Ignore ligands that are not double bonded
                            if any([abs(2 - order) < 1e-7 for order in bond.order]):
                                for ligAtomType in ligand.atomtype:
                                    if ligand.atomtype[0].is_specific_case_of(ATOMTYPES['O']):
                                        correct_atom_list.append('CO')
                                    elif ligand.atomtype[0].is_specific_case_of(ATOMTYPES['S']):
                                        correct_atom_list.append('CS')

                    # remove duplicates from correctAtom:
                    correct_atom_list = list(set(correct_atom_list))
                    for correctAtom in correct_atom_list:
                        tst.append((ATOMTYPES[correctAtom] in atom.atomtype, """
In family {0}, node {1} is missing the atomtype {2} in atom {3} and may be misusing the atomtype Cd, CO, CS, or Cdd.
The following adjList may have atoms in a different ordering than the input file:
{4}""".format(family_name, entry, correctAtom, index + 1, entry.item.to_adjacency_list())))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_unimolecular_groups(self, family_name):
        """
        This test goes through all unimolecular groups that have more than one top level, top level groups
        that overlap with family.reactant are assumed to be backbones(contains the whole reactant molecule)
        and the other top levels are assumedto be endgroups

        the following are format requirements are checked:
        1)endgroup entries hav exactly the same labels as their top level entry
        2)backbone groups have all labels that endgroups have
        3)backbone groups have labels tracing between the endgroups that follow the shortest path
        4)The end subgraph inside each backbone is exactly the same as the top level of the correspodning end tree
        """

        def get_end_from_backbone(_backbone, _end_labels):
            """
            :param _backbone: :class: Entry for a backbone of molecule
            :param _end_labels: Labels in the end groups
            :return: A subgraph representing the end group of the molecule
            """
            # make copy for manipulation
            copy_group = _backbone.item.copy(True)

            # Find the end group atoms
            for _atom in copy_group.atoms:
                if _atom.label in _end_labels:
                    mid_atom = _atom
                    break

            # find the bonds to break
            bonds_to_break = []
            for atom2, bond in mid_atom.bonds.items():
                if atom2.label is None or atom2.label not in _end_labels:
                    bonds_to_break.append(bond)

            for bond in bonds_to_break:
                copy_group.remove_bond(bond)

            # split group into end and backbone fragment
            groups = copy_group.split()

            # verify group was split correctly and identify the correct end group
            _end_labels = set(_end_labels)
            for _group in groups:
                group_labels = set(_atom.label for _atom in _group.atoms)
                group_labels.discard('')
                if _end_labels == group_labels:
                    break
            else:
                raise Exception("Group {0} not split correctly".format(_backbone.label))

            return _group

        #################################################################################
        family = self.database.kinetics.families[family_name]

        backbone = family.get_backbone_roots()[0]

        end_groups = family.get_end_roots()

        end_labels = {}
        for end_group in end_groups:
            labels = []
            for atom in end_group.item.atoms:
                if atom.label:
                    labels.append(atom.label)
            end_labels[end_group] = set(labels)

        # get boundary atoms to test that backbones have labels between end groups
        nose.tools.assert_is_not_none(family.boundary_atoms)

        # set of all end_labels should be backbone label
        backbone_label = set([])
        for end_label in end_labels.values():
            for label in end_label:
                backbone_label.add(label)

        # define types of errors
        a = []  # end groups have too many labels
        b = []  # end group lacks necessary label
        c = []  # backbone missing end group labels
        d = []  # backbone missing labels in between groups
        e = []  # backbone tries to define atoms inside end groups
        for group_name, entry in family.groups.entries.items():
            if isinstance(entry.item, Group):
                group = entry.item
                if backbone in family.ancestors(entry):
                    for atom in group.atoms:
                        if atom.label:
                            present_labels.add(atom.label)
                    # Check C
                    for end_group, labels in end_labels.items():
                        if not labels.issubset(present_labels):
                            c.append([end_group, entry])
                    # check D
                    mid_atoms = [group.get_labeled_atoms(x)[0] for x in family.boundary_atoms]
                    path_atoms = find_shortest_path(mid_atoms[0], mid_atoms[1])
                    for atom in path_atoms:
                        if not atom.label:
                            d.append([backbone, entry])
                            break
                    # check E
                    for end_group, labels in end_labels.items():
                        end_from_backbone = get_end_from_backbone(entry, labels)
                        present_labels = end_from_backbone.get_all_labeled_atoms()
                        present_labels = set(present_labels.keys())
                        if labels == present_labels:
                            if not end_group.item.is_identical(end_from_backbone):
                                e.append([end_group, entry])
                        else:
                            raise Exception(
                                "Group {0} has split into end group {1}, but does not match any root".format(
                                    entry.label, end_from_backbone.to_adjacency_list()))

                else:
                    present_labels = set([])
                    for endNode, labelledAtoms in end_labels.items():
                        if endNode in family.ancestors(entry):
                            for atom in group.atoms:
                                if atom.label:
                                    present_labels.add(atom.label)
                            # Check A
                            if not present_labels.issubset(labelledAtoms):
                                a.append([endNode, entry])
                            # Check B
                            if not labelledAtoms.issubset(present_labels):
                                b.append([endNode, entry])

        # print outputs
        tst = []
        if a:
            s = "These end groups have extra labels that their top level end group do not have:"
            s += "\n [root group, error group]"
            for x in a:
                s += '\n' + str(x)
            tst.append((False, s))
        if b:
            s = "These end groups are missing labels that their top level end group have:"
            s += "\n [root group, error group]"
            for x in b:
                s += '\n' + str(x)
            tst.append((False, s))
        if c:
            s = "These backbone groups are missing labels that are in the end groups:"
            s += "\n [root group, error group]"
            for x in c:
                s += '\n' + str(x)
            tst.append((False, s))
        if d:
            s = "These backbone groups are missing labels along the path atoms:"
            s += "\n [root group, error group]"
            for x in d:
                s += '\n' + str(x)
            tst.append((False, s))
        if e:
            s = "These backbone have end subgraphs that don't match a root:"
            s += "\n [root group, error group]"
            for x in e:
                s += '\n' + str(x)
            tst.append((False, s))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_check_sample_descends_to_group(self, family_name):
        """
        This test first creates a sample :class:Molecule from a :class:Group. Then it checks
        that this molecule hits the original group or a child when it descends down the tree.
        """
        family = self.database.kinetics.families[family_name]

        tst1 = []
        tst2 = []
        tst3 = []
        # ignore any products
        ignore = []
        if not family.own_reverse:
            for product in family.forward_template.products:
                ignore.append(product)
                ignore.extend(product.children)
        else:
            ignore = []

        # If family is backbone archetype, then we need to merge groups before descending
        roots = family.groups.top
        if len(roots) > len(family.forward_template.reactants):
            backbone_roots = family.get_backbone_roots()
            all_backbone_groups = []
            for backboneRoot in backbone_roots:
                all_backbone_groups.extend(family.get_top_level_groups(backboneRoot))
            # list of numbered of labelled atoms for all_backbone_groups
            backbone_sizes = [len(backbone.item.get_all_labeled_atoms()) for backbone in all_backbone_groups]

            # pick a backbone that is two labelled atoms larger than the smallest
            if min(backbone_sizes) + 2 in backbone_sizes:
                backbone_sample = all_backbone_groups[backbone_sizes.index(min(backbone_sizes) + 2)]
            # or if it doesn't exist, pick the largest backbone
            else:
                backbone_sample = all_backbone_groups[backbone_sizes.index(max(backbone_sizes))]
            merges_necessary = True
        else:
            merges_necessary = False

        # If atom has too many benzene rings, we currently have trouble making sample atoms
        skipped = []
        for entryName, entry in family.groups.entries.items():
            if entry in ignore:
                continue
            elif isinstance(entry.item, Group):
                ancestors = family.ancestors(entry)
                if ancestors:
                    root = ancestors[-1]  # top level root will be last one in ancestors
                else:
                    root = entry
                try:
                    if merges_necessary and root not in backbone_roots:  # we may need to merge
                        merged_group = backbone_sample.item.merge_groups(entry.item)
                        sample_molecule = merged_group.make_sample_molecule()
                    else:
                        sample_molecule = entry.item.make_sample_molecule()

                    # test accessibility here
                    atoms = sample_molecule.get_all_labeled_atoms()
                    match = family.groups.descend_tree(sample_molecule, atoms, strict=True, root=root)
                    tst1.append((match, "Group {0} does not match its root node, {1}".format(entryName, root.label)))
                    if tst1[-1][0] is not None:
                        if merges_necessary and root not in backbone_roots:
                            backbone_msg = "\n\nBackbone Group Adjlist:\n" + backbone_sample.label + '\n'
                            backbone_msg += backbone_sample.item.to_adjacency_list()
                        else:
                            backbone_msg = ''
                        tst2.append((entry, [match] + family.groups.ancestors(match), """
In group {0}, a sample molecule made from node {1} returns node {2} when descending the tree.
Sample molecule AdjList:
{3}

Origin Group AdjList:
{4}{5}

Matched group AdjList:
{6}""".format(family_name, entry.label, match.label, sample_molecule.to_adjacency_list(), entry.item.to_adjacency_list(),
              backbone_msg, match.item.to_adjacency_list())))

                except UnexpectedChargeError as e:
                    if merges_necessary and root not in backbone_roots:
                        backbone_msg = "\n\nBackbone Group Adjlist:\n" + backbone_sample.label + '\n'
                        backbone_msg += backbone_sample.item.to_adjacency_list()
                    else:
                        backbone_msg = ''
                    tst3.append((False, """
In family {0}, a sample molecule made from node {1} returns an unexpectedly charged molecule:
Sample molecule AdjList:
{2}

Origin Group AdjList:
{3}{4}""".format(family_name, entry.label, e.graph.to_adjacency_list(), entry.item.to_adjacency_list(), backbone_msg)))

                except ImplicitBenzeneError:
                    skipped.append(entryName)

        # print out entries skipped from exception we can't currently handle
        if skipped:
            print("These entries were skipped because too big benzene rings:")
            for entryName in skipped:
                print(entryName)

        boo = False
        for i in range(len(tst1)):
            if tst1[i][0] is None:
                logging.error(tst1[i][1])
                boo = True
        for i in range(len(tst2)):
            if tst2[i][0] not in tst2[i][1]:
                logging.error(tst2[i][2])
                boo = True
        for i in range(len(tst3)):
            if not tst3[i][0]:
                logging.error(tst3[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def kinetics_check_sample_can_react(self, family_name):
        """
        This test first creates a/some sample :class:Molecule(s) from a/some :class:Group(s). 
        Then it checks that this/these molecule can react according to the recipe.

        It doesn't do every combination of test molecules, but tries each test molecule 
        at least once. Also note that the test molecules are often too small to detect
        some problems, always being the simplest possible instance of a group.
        So unfortunately, problems can still slip through this test.
        """
        family = self.database.kinetics.families[family_name]
        test1 = []

        # ignore any products
        ignore = []
        if not family.own_reverse:
            for product in family.forward_template.products:
                ignore.append(product)
                ignore.extend(product.children)
        else:
            ignore = []

        # If family is backbone archetype, then we need to merge groups before descending
        roots = family.groups.top
        if len(roots) > len(family.forward_template.reactants):
            backbone_roots = family.get_backbone_roots()
            all_backbone_groups = []
            for backboneRoot in backbone_roots:
                all_backbone_groups.extend(family.get_top_level_groups(backboneRoot))
            # list of numbered of labelled atoms for all_backbone_groups
            backbone_sizes = [len(backbone.item.get_all_labeled_atoms()) for backbone in all_backbone_groups]

            # pick a backbone that is two labelled atoms larger than the smallest
            if min(backbone_sizes) + 2 in backbone_sizes:
                backbone_sample = all_backbone_groups[backbone_sizes.index(min(backbone_sizes) + 2)]
            # or if it doesn't exist, pick the largest backbone
            else:
                backbone_sample = all_backbone_groups[backbone_sizes.index(max(backbone_sizes))]
            merges_necessary = True
        else:
            merges_necessary = False

        # If atom has too many benzene rings, we currently have trouble making sample atoms
        skipped = []
        sample_reactants = defaultdict(list) # the keys will be the root nodes, the values a list of samples
        for entryName, entry in family.groups.entries.items():
            if entry in ignore:
                continue
            elif isinstance(entry.item, Group):
                ancestors = family.ancestors(entry)
                if ancestors:
                    root = ancestors[-1]  # top level root will be last one in ancestors
                else:
                    root = entry
                try:
                    if merges_necessary and root not in backbone_roots:  # we may need to merge
                        merged_group = backbone_sample.item.merge_groups(entry.item)
                        sample_molecule = merged_group.make_sample_molecule()
                        # store the sample molecule for later testing
                        if not family.is_molecule_forbidden(sample_molecule):
                            sample_reactants[backboneRoot].append(sample_molecule)
                    else:
                        sample_molecule = entry.item.make_sample_molecule()
                        # store the sample molecule for later testing
                        if not family.is_molecule_forbidden(sample_molecule):
                            sample_reactants[root].append(sample_molecule)
                except UnexpectedChargeError as e:
                    if merges_necessary and root not in backbone_roots:
                        backbone_msg = "\n\nBackbone Group Adjlist:\n" + backbone_sample.label + '\n'
                        backbone_msg += backbone_sample.item.to_adjacency_list()
                    else:
                        backbone_msg = ''
                    tst3.append((False, """
In family {0}, a sample molecule made from node {1} returns an unexpectedly charged molecule:
Sample molecule AdjList:
{2}

Origin Group AdjList:
{3}{4}""".format(family_name, entry.label, e.graph.to_adjacency_list(), entry.item.to_adjacency_list(), backbone_msg)))

                except ImplicitBenzeneError:
                    skipped.append(entryName)
        def make_error_message(reactants, message=''):
            "Give it the list of reactant Molecules and an optional message."
            output = f"Error in family {family_name} when reacting "
            output += ' + '.join(s.to_smiles() for s in reactants)
            output += f". {message}\n"
            for s in reactants:
                output += "\n" + s.to_adjacency_list(label=s.to_smiles())
            return output

        if len(sample_reactants) == 1 == len(family.forward_template.reactants):
            reactants = list(sample_reactants.values())[0]
            for reactant in reactants:
                try:
                    products = family.apply_recipe([reactant])
                except Exception as e:
                    test1.append(make_error_message([reactant],
                          message=f"During apply_recipe had an {type(e)!s}: {e!s}"))
                    continue
                if products is None:
                    test1.append(make_error_message([reactant],
                        message="apply_recipe returned None, indicating wrong number of products or a charged product."))
                    continue
                for molecule in products:
                    # Just check none of this throws errors
                    species = rmgpy.species.Species(index=1,molecule=[molecule])
                    species.generate_resonance_structures()
        elif len(sample_reactants) == 1 and len(family.forward_template.reactants) == 2:
            # eg. Bimolec_Hydroperoxide_Decomposition and Peroxyl_Disproportionation
            # use the same group twice.
            reactant_lists = [sample_reactants[k] for k in family.forward_template.reactants ]
            pairs = zip(*reactant_lists)
            for reactant1, reactant2 in pairs:
                try:
                    products = family.apply_recipe([reactant1, reactant2])
                except Exception as e:
                    test1.append(make_error_message([reactant1, reactant2],
                          message=f"During apply_recipe had an {type(e)!s}: {e!s}"))
                    continue
                if products is None:
                    test1.append(make_error_message([reactant1, reactant2],
                        message="apply_recipe returned None, indicating wrong number of products or a charged product."))
                    continue
                for molecule in products:
                    # Just check none of this throws errors
                    species = rmgpy.species.Species(index=1,molecule=[molecule])
                    species.generate_resonance_structures()
        elif len(sample_reactants) == 2:
            sr = list(sample_reactants.values())
            # Every combination may be prohibitively slow (N*M reactions),
            # so we loop through to ensure at least each node is used once (max(N,M) reactions)
            if len(sr[0]) > len(sr[1]):
                pairs = zip(sr[0], itertools.cycle(sr[1]))
            else:
                pairs = zip(itertools.cycle(sr[0]), sr[1])
            for reactant1, reactant2 in pairs:
                try:
                    products = family.apply_recipe([reactant1, reactant2])
                except Exception as e:
                    test1.append(make_error_message([reactant1, reactant2],
                          message=f"During apply_recipe had an {type(e)!s}: {e!s}"))
                    continue
                if products is None:
                    test1.append(make_error_message([reactant1, reactant2],
                        message="apply_recipe returned None, indicating wrong number of products or a charged product."))
                    continue
                for molecule in products:
                    # Just check none of this throws errors
                    species = rmgpy.species.Species(index=1,molecule=[molecule])
                    species.generate_resonance_structures()
        elif len(sample_reactants) == 3:
            sr = list(sample_reactants.values())
            # Every combination may be prohibitively slow (N*M*L reactions),
            # so we loop through to ensure at least each node is used once (max(N,M,L) reactions)
            longest = np.argmax(map(len,sr))
            if longest == 0:
                triplets = zip(sr[0], itertools.cycle(sr[1]), itertools.cycle(sr[2]))
            elif longest == 1:
                triplets = zip(itertools.cycle(sr[0]), sr[1], itertools.cycle(sr[2]))
            else:
                assert longest == 2
                triplets = zip(itertools.cycle(sr[0]), itertools.cycle(sr[1]), sr[2])
            for reactant1, reactant2, reactant3 in triplets:
                try:
                    products = family.apply_recipe([reactant1, reactant2, reactant3])
                except Exception as e:
                    test1.append(make_error_message([reactant1, reactant2, reactant3],
                          message=f"During apply_recipe had an {type(e)!s}: {e!s}"))
                    continue
                if products is None:
                    test1.append(make_error_message([reactant1, reactant2, reactant3],
                        message="apply_recipe returned None, indicating wrong number of products or a charged product."))
                    continue
                for molecule in products:
                    # Just check none of this throws errors
                    species = rmgpy.species.Species(index=1,molecule=[molecule])
                    species.generate_resonance_structures()
        else:
            raise ValueError(f"Family had {len(sample_reactants)} reactants?: "
                             f"{', '.join(map(str,sample_reactants.keys())) }")

        # print out entries skipped from exception we can't currently handle
        if skipped:
            print("These entries were skipped because of an ImplicitBenzeneError:")
            for entryName in skipped:
                print(entryName)

        boo = False
        for err in test1:
            logging.error(err)
            boo = True
        if boo:
            raise ValueError("Error Occurred. See log for details.")
    def check_surface_thermo_groups_have_surface_attributes(self, group_name, group):
        """
        Tests that each entry in the surface thermo groups has a 'metal' and 'facet' attribute, 
        describing which metal the data came from.
        """
        failed = False
        for entry in group.entries.values():
            if isinstance(entry.data, rmgpy.thermo.thermodata.ThermoData):
                if 'Pt' in group_name:
                    if entry.metal is not 'Pt':
                        logging.error(f'Expected {entry} metal attribute in {group_name} group to match Pt, but was {entry.metal}')
                        failed = True
                if '111' in group_name:
                    if entry.facet is not '111':
                        logging.error(f'Expected {entry} facet attribute in {group_name} group to match 111, but was {entry.facet}')
                        failed = True
                if not entry.metal:
                    logging.error(f'Expected a metal attribute for {entry} in {group_name} group but found {entry.metal!r}')
                    failed = True
                if not entry.facet:
                    logging.error(f'Expected a facet attribute for {entry} in {group_name} group but found {entry.facet!r}')
                    failed = True
        if failed:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def check_surface_thermo_libraries_have_surface_attributes(self, library_name, library):
        """
        Test that each entry in the surface thermo database has a 'metal' and 'facet' attribute,
        describing which metal the data came from.
        """
        failed = False
        for entry in library.entries.values():
            if 'Pt' in library_name:
                if entry.metal is not 'Pt':
                    logging.error(f'Expected {entry} metal attribute in {library_name} library to match Pt, but was {entry.metal}')
                    failed = True
            if 'Ni' in library_name:
                if entry.metal is not 'Ni':
                    logging.error(f'Expected {entry} metal attribute in {library_name} library to match Ni, but was {entry.metal}')
                    failed = True
            if '111' in library_name:
                if entry.facet is not '111':
                    logging.error(f'Expected {entry} facet attribute in {library_name} library to match 111, but was {entry.facet}')
                    failed = True
            if not entry.metal:
                logging.error(f'Expected a metal attribute for {entry} in {library} library but found {entry.metal!r}')
                failed = True
            if not entry.facet:
                logging.error(f'Expected a facet attribute for {entry} in {library} library but found {entry.facet!r}')
                failed = True
            if failed:
                raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def general_check_nodes_found_in_tree(self, group_name, group):
        """
        This test checks whether nodes are found in the tree, with proper parents.
        """
        for node_name, node_group in group.entries.items():
            ascend_parent = node_group
            # Check whether the node has proper parents unless it is the top reactant or product node
            tst1 = []
            tst2 = []
            tst3 = []
            while ascend_parent not in group.top:
                child = ascend_parent
                ascend_parent = ascend_parent.parent
                tst1.append((ascend_parent is not None,
                             "Node {node} in {group} group was found in the tree without a proper parent.".format(
                                 node=child, group=group_name)))
                if tst1[-1] is not None:
                    tst2.append((child in ascend_parent.children,
                                 "Node {node} in {group} group was found in the tree without a proper parent.".format(
                                     node=node_name, group=group_name)))
                    tst3.append((child is ascend_parent,
                                 "Node {node} in {group} is a parent to itself".format(
                                     node=node_name, group=group_name)))

        boo = False
        for i in range(len(tst1)):
            if not tst1[i][0]:
                logging.error(tst1[i][1])
                boo = True
        for i in range(len(tst2)):
            if not tst2[i][0]:
                logging.error(tst2[i][1])
                boo = True
            if tst3[i][0]:
                logging.error(tst3[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_check_groups_nonidentical(self, group_name, group):
        """
        This test checks whether nodes found in the group are nonidentical.
        """
        entries_copy = copy(group.entries)
        tst = []
        for node_name, node_group in group.entries.items():
            del entries_copy[node_name]
            for node_name_other, node_group_other in entries_copy.items():
                group.match_node_to_node(node_group, node_group_other)
                tst.append((group.match_node_to_node(node_group, node_group_other),
                            "Node {node} in {group} group was found to be identical to node {node_other}".format(
                                node=node_name, group=group_name, node_other=node_name_other)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_check_child_parent_relationships(self, group_name, group):
        """
        This test checks that nodes' parent-child relationships are correct in the database.
        """
        tst1 = []
        tst2 = []
        for node_name, child_node in group.entries.items():
            # top nodes and product nodes don't have parents by definition, so they get an automatic pass:
            if child_node in group.top:
                continue
            parent_node = child_node.parent
            # Check whether the node has proper parents unless it is the top reactant or product node
            # The parent should be more general than the child
            tst1.append((group.match_node_to_child(parent_node, child_node),
                         "In {group} group, node {parent} is not a proper parent of its child {child}.".format(
                             group=group_name, parent=parent_node, child=node_name)))

            # check that parentNodes which are LogicOr do not have an ancestor that is a Group
            # If it does, then the child_node must also be a child of the ancestor
            if isinstance(parent_node.item, LogicOr):
                ancestor_node = parent_node
                while ancestor_node not in group.top and isinstance(ancestor_node.item, LogicOr):
                    ancestor_node = ancestor_node.parent
                if isinstance(ancestor_node.item, Group) and tst1[-1][0]:
                    tst2.append((group.match_node_to_child(ancestor_node, child_node),
                                 "In {group} group, node {ancestor} is not a proper ancestor of its child {child}."
                                 "".format(group=group_name, ancestor=ancestor_node, child=node_name)))

        boo = False
        for i in range(len(tst1)):
            if not tst1[i][0]:
                logging.error(tst1[i][1])
                boo = True
        for i in range(len(tst2)):
            if not tst2[i][0]:
                logging.error(tst2[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_check_siblings_for_parents(self, group_name, group):
        """
        This test checks that siblings in a tree are not actually parent/child.

        For example in a tree:

        L1. A
            L2. B
            L2. C

        This tests that C is not a child of B, which would make C inaccessible because
        we always match B first.

        We do not check that B is not a child of C becausethat does not cause accessibility
        problems and may actually be necessary in some trees. For example, in the polycyclic
        thermo groups B might be a tricyclic and C a bicyclic parent. Currently there is no
        way to writes a bicyclic group that excludes an analogous tricyclic.
        """
        tst = []
        for nodeName, node in group.entries.items():
            for index, child1 in enumerate(node.children):
                for child2 in node.children[index + 1:]:
                    tst.append((group.match_node_to_child(child1, child2),
                                "In {0} group, node {1} is a parent of {2}, but they are written as siblings.".format(
                                    group_name, child1, child2)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_check_cd_atom_type(self, group_name, group):
        """
        This test checks that groups containing Cd, CO, CS and Cdd atomtypes are used
        correctly according to their strict definitions
        """
        target_label = ['Cd', 'CO', 'CS', 'Cdd']
        target_atom_types = [ATOMTYPES[x] for x in target_label]
        tst = []
        for entry_name, entry in group.entries.items():
            if isinstance(entry.item, Group):
                for index, atom in enumerate(entry.item.atoms):
                    for atomtype1 in atom.atomtype:
                        if atomtype1 in target_atom_types:
                            break
                    else:
                        # If Cd not found in atomTypes, go to next atom
                        continue
                    # figure out what the correct atomtype is
                    correct_atom_list = []
                    num_of_d_bonds = sum(
                        [1 if x.order[0] is 'D' and len(x.order) == 1 else 0 for x in atom.bonds.values()])
                    if num_of_d_bonds == 2:
                        correct_atom_list.append('Cdd')
                    elif num_of_d_bonds == 1:
                        for ligand, bond in atom.bonds.items():
                            # Ignore ligands that are not double bonded
                            if any([abs(2 - order) < 1e-7 for order in bond.order]):
                                for lig_atom_type in ligand.atomtype:
                                    if ligand.atomtype[0].is_specific_case_of(ATOMTYPES['O']):
                                        correct_atom_list.append('CO')
                                    elif ligand.atomtype[0].is_specific_case_of(ATOMTYPES['S']):
                                        correct_atom_list.append('CS')

                    # remove duplicates from correctAtom:
                    correct_atom_list = list(set(correct_atom_list))
                    for correctAtom in correct_atom_list:
                        tst.append((ATOMTYPES[correctAtom] in atom.atomtype, """
In group {0}, node {1} is missing the atomtype {2} in atom {3} and may be misusing the atomtype Cd, CO, CS, or Cdd.
The following adjList may have atoms in a different ordering than the input file:
{4}""".format(group_name, entry, correctAtom, index + 1, entry.item.to_adjacency_list())))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_check_sample_descends_to_group(self, group_name, group):
        """
        This test first creates a sample :class:Molecule from a :class:Group. Then it checks
        that this molecule hits the original group or a child when it descends down the tree.
        """

        skipped = []
        tst1 = []
        tst2 = []
        tst3 = []
        for entryName, entry in group.entries.items():
            try:
                if isinstance(entry.item, Group):
                    try:
                        sample_molecule = entry.item.make_sample_molecule()
                    except:
                        logging.error("Problem making sample molecule for group {}\n{}".format(
                            entryName, entry.item.to_adjacency_list()))
                        raise

                    atoms = sample_molecule.get_all_labeled_atoms()
                    match = group.descend_tree(sample_molecule, atoms, strict=True)
                    tst1.append((match, "Group {0} does not match its root node, {1}".format(entryName, group.top[0])))
                    tst2.append((entry, [match] + group.ancestors(match), """
In group {0}, a sample molecule made from node {1} returns node {2} when descending the tree.
Sample molecule AdjList:
{3}

Origin Group AdjList:
{4}

Matched group AdjList:
{5}
""".format(group_name,
           entry,
           match,
           sample_molecule.to_adjacency_list(),
           entry.item.to_adjacency_list(),
           match.item.to_adjacency_list())))

            except UnexpectedChargeError as e:
                tst3.append((False, """
In family {0}, a sample molecule made from node {1} returns an unexpectedly charged molecule:
Sample molecule AdjList:
{2}

Origin Group AdjList:
{3}""".format(group_name, entry.label, e.graph.to_adjacency_list(), entry.item.to_adjacency_list())))

            except ImplicitBenzeneError:
                skipped.append(entryName)

        # print out entries skipped from exception we can't currently handle
        if skipped:
            print("These entries were skipped because too big benzene rings:")
            for entryName in skipped:
                print(entryName)

        boo = False
        for i in range(len(tst1)):
            if tst1[i][0] is None:
                logging.error(tst1[i][1])
                boo = True
            if tst2[i][0] not in tst2[i][1]:
                logging.error(tst2[i][2])
                boo = True
        for i in range(len(tst3)):
            if not tst3[i][0]:
                logging.error(tst3[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")


if __name__ == '__main__':
    nose.run(argv=[__file__, '-v', '--nologcapture'], defaultTest=__name__)
