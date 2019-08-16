"""
This scripts runs tests on the database
"""
from __future__ import division, print_function

import numpy as np
import logging

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from copy import copy
import rmgpy.kinetics
import quantities as pq
from rmgpy.data.base import LogicOr
from rmgpy.molecule import Group, ImplicitBenzeneError, UnexpectedChargeError
from rmgpy.molecule.atomtype import atomTypes
from rmgpy.molecule.pathfinder import find_shortest_path

import nose
import nose.tools


class TestDatabase(object):  # cannot inherit from unittest.TestCase if we want to use nose test generators
    """
    Contains unit tests for the database for rigorous error checking.
    """
    @classmethod
    def setUpClass(cls):
        """
        Load the database before running the tests.
        """
        databaseDirectory = settings['database.directory']
        cls.database = RMGDatabase()
        cls.database.load(databaseDirectory, kineticsFamilies='all')

    # These are generators, that call the methods below.
    def test_kinetics(self):
        for family_name, family in self.database.kinetics.families.items():

            test = lambda x: self.kinetics_checkCorrectNumberofNodesInRules(family_name)
            test_name = "Kinetics family {0}: rules have correct number of nodes?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, None

            test = lambda x: self.kinetics_checkNodesInRulesFoundInGroups(family_name)
            test_name = "Kinetics family {0}: rules' nodes exist in the groups?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, None

            test = lambda x: self.kinetics_checkGroupsFoundInTree(family_name)
            test_name = "Kinetics family {0}: groups are in the tree with proper parents?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, None

            test = lambda x: self.kinetics_checkGroupsNonidentical(family_name)
            test_name = "Kinetics family {0}: groups are not identical?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_checkChildParentRelationships(family_name)
            test_name = "Kinetics family {0}: parent-child relationships are correct?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_checkSiblingsForParents(family_name)
            test_name = "Kinetics family {0}: sibling relationships are correct?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_checkCdAtomType(family_name)
            test_name = "Kinetics family {0}: Cd, CS, CO, and Cdd atomtype used correctly?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            test = lambda x: self.kinetics_checkReactantAndProductTemplate(family_name)
            test_name = "Kinetics family {0}: reactant and product templates correctly defined?".format(family_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, family_name

            #these families have some sort of difficulty which prevents us from testing accessibility right now
            difficultFamilies = ['Diels_alder_addition', 'Intra_R_Add_Exocyclic', 'Intra_R_Add_Endocyclic']
            generatedTrees = ["R_Recombination"]

            if len(family.forwardTemplate.reactants) < len(family.groups.top) and family_name not in difficultFamilies:
                test = lambda x: self.kinetics_checkUnimolecularGroups(family_name)
                test_name = "Kinetics family {0} check that unimolecular group is formatted correctly?".format(family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

            if family_name not in difficultFamilies and family_name not in generatedTrees:
                test = lambda x: self.kinetics_checkSampleDescendsToGroup(family_name)
                test_name = "Kinetics family {0}: Entry is accessible?".format(family_name)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, family_name

            for depository in family.depositories:

                test = lambda x: self.kinetics_checkAdjlistsNonidentical(depository)
                test_name = "Kinetics depository {0}: check adjacency lists are nonidentical?".format(depository.label)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, depository.label

                test = lambda x: self.kinetics_checkRateUnitsAreCorrect(depository, tag='depository')
                test_name = "Kinetics depository {0}: check rates have correct units?".format(depository.label)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, depository.label

        for library_name, library in self.database.kinetics.libraries.items():

            test = lambda x: self.kinetics_checkAdjlistsNonidentical(library)
            test_name = "Kinetics library {0}: check adjacency lists are nonidentical?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

            test = lambda x: self.kinetics_checkRateUnitsAreCorrect(library)
            test_name = "Kinetics library {0}: check rates have correct units?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

            test = lambda x: self.kinetics_checkLibraryRatesAreReasonable(library)
            test_name = "Kinetics library {0}: check rates are reasonable?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name

    def test_thermo(self):
        for group_name, group in self.database.thermo.groups.items():
            test = lambda x: self.general_checkNodesFoundInTree(group_name, group)
            test_name = "Thermo groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkGroupsNonidentical(group_name, group)
            test_name = "Thermo groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkChildParentRelationships(group_name, group)
            test_name = "Thermo groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSiblingsForParents(group_name, group)
            test_name = "Thermo groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkCdAtomType(group_name, group)
            test_name = "Thermo groups {0}: Cd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSampleDescendsToGroup(group_name, group)
            test_name = "Thermo groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

    def test_solvation(self):
        for group_name, group in self.database.solvation.groups.items():
            test = lambda x: self.general_checkNodesFoundInTree(group_name, group)
            test_name = "Solvation groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkGroupsNonidentical(group_name, group)
            test_name = "Solvation groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkChildParentRelationships(group_name, group)
            test_name = "Solvation groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSiblingsForParents(group_name, group)
            test_name = "Solvation groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkCdAtomType(group_name, group)
            test_name = "Solvation groups {0}: Cd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSampleDescendsToGroup(group_name, group)
            test_name = "Solvation groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

    def test_statmech(self):
        for group_name, group in self.database.statmech.groups.items():
            test = lambda x: self.general_checkNodesFoundInTree(group_name, group)
            test_name = "Statmech groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkGroupsNonidentical(group_name, group)
            test_name = "Statmech groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkChildParentRelationships(group_name, group)
            test_name = "Statmech groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSiblingsForParents(group_name, group)
            test_name = "Statmech groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkCdAtomType(group_name, group)
            test_name = "Statmech groups {0}: Cd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSampleDescendsToGroup(group_name, group)
            test_name = "Statmech groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

    def test_transport(self):
        for group_name, group in self.database.transport.groups.items():
            test = lambda x: self.general_checkNodesFoundInTree(group_name, group)
            test_name = "Transport groups {0}: nodes are in the tree with proper parents?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkGroupsNonidentical(group_name, group)
            test_name = "Transport groups {0}: nodes are nonidentical?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkChildParentRelationships(group_name, group)
            test_name = "Transport groups {0}: parent-child relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSiblingsForParents(group_name, group)
            test_name = "Transport groups {0}: sibling relationships are correct?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkCdAtomType(group_name, group)
            test_name = "Transport groups {0}: Cd, CS, CO, and Cdd atomtype used correctly?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name

            test = lambda x: self.general_checkSampleDescendsToGroup(group_name, group)
            test_name = "Transport groups {0}: Entry is accessible?".format(group_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, group_name


    # These are the actual tests, that don't start with a "test_" name:
    def kinetics_checkCorrectNumberofNodesInRules(self, family_name):
        """
        This test ensures that each rate rule contains the proper number of nodes according to the family it originates.
        """
        family = self.database.kinetics.families[family_name]
        expectedNumberNodes = len(family.getRootTemplate())
        tst = []
        for label, entries in family.rules.entries.items():
            for entry in entries:
                nodes = label.split(';')
                tst.append((len(nodes), expectedNumberNodes, "Wrong number of groups or semicolons in family {family} rule {entry}.  Should be {num_nodes}".format(family=family_name, entry=entry, num_nodes=expectedNumberNodes)))

        boo = False
        for item in tst:
            if item[0] != item[1]:
                boo = True
                logging.error(item[2])
        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkNodesInRulesFoundInGroups(self, family_name):
        """
        This test ensures that each rate rule contains nodes that exist in the groups and that they match the order of the forwardTemplate.
        """
        family = self.database.kinetics.families[family_name]

        # List of the each top node's descendants (including the top node)
        topDescendants = []
        for topNode in family.getRootTemplate():
            nodes = [topNode]
            nodes.extend(family.groups.descendants(topNode))
            topDescendants.append(nodes)

        topGroupOrder = ';'.join(topNode.label for topNode in family.getRootTemplate())
        tst1 = []
        tst2 = []
        for label, entries in family.rules.entries.items():
            for entry in entries:
                nodes = label.split(';')
                for i, node in enumerate(nodes):
                    tst1.append((node in family.groups.entries, "In {family} family, no group definition found for label {label} in rule {entry}".format(family=family_name, label=node, entry=entry)))
                    tst2.append((family.groups.entries[node] in topDescendants[i], "In {family} family, rule {entry} was found with groups out of order.  The correct order for a rule should be subgroups of {top}.".format(family=family_name, entry=entry, top=topGroupOrder)))
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

    def kinetics_checkGroupsFoundInTree(self, family_name):
        """
        This test checks whether groups are found in the tree, with proper parents.
        """
        family = self.database.kinetics.families[family_name]
        tst = []
        tst1 = []
        tst2 = []
        tst3 = []
        for nodeName, nodeGroup in family.groups.entries.items():
            tst.append(('[' in nodeName or ']' in nodeName, "Group {group} in {family} family contains square brackets [ ] in the label, which are not allowed.".format(group=nodeName, family=family_name)))
            ascendParent = nodeGroup

            # Check whether the node has proper parents unless it is the top reactant or product node
            while ascendParent not in family.groups.top and ascendParent not in family.forwardTemplate.products:
                child = ascendParent
                ascendParent = ascendParent.parent
                tst1.append((ascendParent is not None, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=child, family=family_name)))
                tst2.append((child in ascendParent.children, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=nodeName, family=family_name)))
                tst3.append((child is ascendParent, "Group {group} in {family} family is a parent to itself".format(group=nodeName, family=family_name)))

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

    def kinetics_checkGroupsNonidentical(self, family_name):
        """
        This test checks that the groups are non-identical.
        """
        from rmgpy.data.base import Database
        originalFamily = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = originalFamily.groups.entries
        entriesCopy = copy(family.entries)
        tst = []
        for nodeName, nodeGroup in family.entries.items():
            del entriesCopy[nodeName]
            for nodeNameOther, nodeGroupOther in entriesCopy.items():
                tst.append((family.matchNodeToNode(nodeGroup, nodeGroupOther), "Group {group} in {family} family was found to be identical to group {groupOther}".format(group=nodeName, family=family_name, groupOther=nodeNameOther)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkChildParentRelationships(self, family_name):
        """
        This test checks that groups' parent-child relationships are correct in the database.
        """
        from rmgpy.data.base import Database
        originalFamily = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = originalFamily.groups.entries
        tst = []
        for nodeName, childNode in family.entries.items():
            #top nodes and product nodes don't have parents by definition, so they get an automatic pass:
            if childNode in originalFamily.groups.top or childNode in originalFamily.forwardTemplate.products: continue
            parentNode = childNode.parent

            if parentNode is None:
                # This is a mistake in the database, but it should be caught by kinetics_checkGroupsFoundInTree
                # so rather than report it twice or crash, we'll just silently carry on to the next node.
                continue
            elif parentNode in originalFamily.forwardTemplate.products:
                #This is a product node made by training reactions which we do not need to heck
                continue
            # Check whether the node has proper parents unless it is the top reactant or product node
            # The parent should be more general than the child
            tst.append((family.matchNodeToChild(parentNode, childNode),
                            "In {family} family, group {parent} is not a proper parent of its child {child}.".format(family=family_name, parent=parentNode, child=nodeName))
)
            #check that parentNodes which are LogicOr do not have an ancestor that is a Group
            #If it does, then the childNode must also be a child of the ancestor
            if isinstance(parentNode.item, LogicOr):
                ancestorNode = parentNode
                while ancestorNode not in originalFamily.groups.top and isinstance(ancestorNode.item, LogicOr):
                    ancestorNode = ancestorNode.parent
                if isinstance(ancestorNode.item, Group):
                    tst.append((family.matchNodeToChild(ancestorNode, childNode),
                                    "In {family} family, group {ancestor} is not a proper ancestor of its child {child}.".format(family=family_name, ancestor=ancestorNode, child=nodeName)))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkSiblingsForParents(self, family_name):
        """
        This test checks that siblings in a tree are not actually parent/child

        See general_checkSiblingsForParents comments for more detailed description
        of the test.
        """
        from rmgpy.data.base import Database
        originalFamily = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = originalFamily.groups.entries
        tst = []
        for nodeName, node in family.entries.items():
            #Some families also construct a 2-level trees for the products
            #(root with all entries down one level) We don't care about this
            #tree as it is not used in searching, so we ignore products
            if node in originalFamily.forwardTemplate.products: continue
            for index, child1 in enumerate(node.children):
                for child2 in node.children[index+1:]:
                    tst.append((family.matchNodeToChild(child1, child2),
                                            "In family {0}, node {1} is a parent of {2}, but they are written as siblings.".format(family_name, child1, child2)))
        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkAdjlistsNonidentical(self, database):
        """
        This test checks whether adjacency lists of reactants in a KineticsDepository or KineticsLibrary database object are nonidentical.
        """
        speciesDict = {}
        entries = database.entries.values()
        for entry in entries:
            for reactant in entry.item.reactants:
                if reactant.label not in speciesDict:
                    speciesDict[reactant.label] = reactant

            for product in entry.item.products:
                if product.label not in speciesDict:
                    speciesDict[product.label] = product

        tst = []
        # Go through all species to make sure they are nonidentical
        speciesList = speciesDict.values()
        labeledAtoms = [species.molecule[0].getLabeledAtoms() for species in speciesList]
        for i in range(len(speciesList)):
            for j in range(i+1,len(speciesList)):
                    initialMap = {}
                    try:
                        atom_labels = set(list(labeledAtoms[i].keys()) +
                                                list(labeledAtoms[j].keys()))
                        for atomLabel in atom_labels:
                            initialMap[labeledAtoms[i][atomLabel]] = labeledAtoms[j][atomLabel]
                    except KeyError:
                        # atom labels did not match, therefore not a match
                        continue
                    tst.append((speciesList[i].molecule[0].isIsomorphic(speciesList[j].molecule[0], initialMap), "Species {0} and species {1} in {2} database were found to be identical.".format(speciesList[i].label,speciesList[j].label,database.label)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkRateUnitsAreCorrect(self, database, tag='library'):
        """
        This test ensures that every reaction has acceptable units on the A factor.
        """
        boo = False

        dimensionalities = {
            1: (1 / pq.s).dimensionality ,
            2: (pq.m**3 / pq.mole / pq.s).dimensionality,
            3: ((pq.m**6) / (pq.mole**2) / pq.s).dimensionality,
        }

        for entry in database.entries.values():
            k = entry.data
            rxn = entry.item
            molecularity = len(rxn.reactants)
            surface_reactants = sum([1 for s in rxn.reactants if s.containsSurfaceSite()])
            try:

                if isinstance(k, rmgpy.kinetics.StickingCoefficient):
                    "Should be dimensionless"
                    A = k.A
                    if A.units:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid units {3}'.format(rxn, tag, database.label, A.units))
                elif isinstance(k, rmgpy.kinetics.SurfaceArrhenius):
                    A = k.A
                    expected = copy(dimensionalities[molecularity])
                    # for each surface reactant but one, switch from (m3/mol) to (m2/mol)
                    expected[pq.m] -= (surface_reactants-1)
                    if pq.Quantity(1.0, A.units).simplified.dimensionality != expected :
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid units {3}'.format(rxn, tag, database.label, A.units))
                elif isinstance(k, rmgpy.kinetics.Arrhenius): # (but not SurfaceArrhenius, which came first)
                    A = k.A
                    if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid units {3}'.format(rxn, tag, database.label, A.units))
                elif isinstance(k, (rmgpy.kinetics.Lindemann, rmgpy.kinetics.Troe )):
                    A = k.arrheniusHigh.A
                    if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid high-pressure limit units {3}'.format(rxn, tag, database.label, A.units))
                elif isinstance(k, (rmgpy.kinetics.Lindemann, rmgpy.kinetics.Troe, rmgpy.kinetics.ThirdBody)):
                    A = k.arrheniusLow.A
                    if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity+1]:
                        boo = True
                        logging.error('Reaction {0} from {1} {2}, has invalid low-pressure limit units {3}'.format(rxn, tag, database.label, A.units))
                elif hasattr(k, 'highPlimit') and k.highPlimit is not None:
                    A = k.highPlimit.A
                    if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity-1]:
                        boo = True
                        logging.error(
                            'Reaction {0} from {1} {2}, has invalid high-pressure limit units {3}'.format(rxn, tag, database.label, A.units))
                elif isinstance(k, rmgpy.kinetics.MultiArrhenius):
                    for num, arrhenius in enumerate(k.arrhenius):
                        A = arrhenius.A
                        if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity]:
                            boo = True
                            logging.error(
                                'Reaction {0} from {1} {2}, has invalid units {3} on rate expression {4}'.format(
                                    rxn, tag, database.label, A.units, num + 1)
                            )

                elif isinstance(k, rmgpy.kinetics.PDepArrhenius):
                    for pa, arrhenius in zip(k.pressures.value_si, k.arrhenius):
                        P = rmgpy.quantity.Pressure(1, k.pressures.units)
                        P.value_si = pa

                        if isinstance(arrhenius, rmgpy.kinetics.MultiArrhenius):
                            # A PDepArrhenius may have MultiArrhenius within it
                            # which is distinct (somehow) from MultiPDepArrhenius
                            for num, arrhenius2 in enumerate(arrhenius.arrhenius):
                                A = arrhenius2.A
                                if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity]:
                                    boo = True
                                    logging.error(
                                        'Reaction {0} from {1} {2}, has invalid units {3} on {4!r} rate expression {5}'.format(
                                            rxn, tag, database.label, A.units, P, num + 1)
                                    )
                        else:
                            A = arrhenius.A
                            if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[molecularity]:
                                boo = True
                                logging.error(
                                    'Reaction {0} from {1} {2}, has invalid {3!r} units {4}'.format(
                                        rxn, tag, database.label, P, A.units)
                                )

                elif isinstance(k, rmgpy.kinetics.MultiPDepArrhenius):
                    for num, k2 in enumerate(k.arrhenius):
                        for pa, arrhenius in zip(k2.pressures.value_si, k2.arrhenius):
                            P = rmgpy.quantity.Pressure(1, k2.pressures.units)
                            P.value_si = pa
                            if isinstance(arrhenius, rmgpy.kinetics.MultiArrhenius):
                                # A MultiPDepArrhenius may have MultiArrhenius within it
                                for arrhenius2 in arrhenius.arrhenius:
                                    A = arrhenius2.A
                                    if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[
                                        molecularity]:
                                        boo = True
                                        logging.error(
                                            'Reaction {0} from {1} {2}, has invalid units {3} on {4!r} rate expression {5!r}'.format(
                                                rxn, tag, database.label, A.units, P, arrhenius2)
                                        )
                            else:
                                A = arrhenius.A
                                if pq.Quantity(1.0, A.units).simplified.dimensionality != dimensionalities[
                                    molecularity]:
                                    boo = True
                                    logging.error(
                                        'Reaction {0} from {1} {2}, has invalid {3!r} units {4} in rate expression {5}'.format(
                                        rxn, tag, database.label, P, A.units, num)
                                    )


                elif isinstance(k, rmgpy.kinetics.Chebyshev):
                    if pq.Quantity(1.0, k.kunits).simplified.dimensionality != dimensionalities[molecularity]:
                        boo = True
                        logging.error(
                            'Reaction {0} from {1} {2}, has invalid units {3}'.format(
                                rxn, tag, database.label, k.kunits)
                        )

                else:
                    logging.warning('Reaction {0} from {1} {2}, did not have units checked.'.format(rxn, tag, database.label))
            except:
                logging.error("Error when checking units on reaction {0} from {1} {2} with rate expression {3!r}.".format(rxn, tag, database.label, k))
                raise
        if boo:
            raise ValueError('{0} {1} has some incorrect units'.format(tag.capitalize(), database.label))

    def kinetics_checkLibraryRatesAreReasonable(self, library):
        """
        This test ensures that every library reaction has reasonable kinetics at 1000 K, 1 bar
        """
        T = 1000.0 #K
        P = 100000.0 # 1 bar in Pa
        Na = 6.022*10**23 #molecules/mol
        mHrad = 1.6737236e-27 #kg
        Hrad_diam = 2.4e-10 #m
        kB = 1.38064852e-23 #m2 * kg *s^-2 * K^-1
        h = 6.62607004e-34 #m2 kg / s
        boo = False
        TST_limit = (kB*T)/h
        collision_limit = Na*np.pi*Hrad_diam**2*np.sqrt(8*kB*T/(np.pi*mHrad/2))
        for entry in library.entries.values():
            if entry.item.isSurfaceReaction():
                # Don't check surface reactions
                continue
            k = entry.data.getRateCoefficient(T,P)
            rxn = entry.item
            if k < 0:
                boo = True
                logging.error('library reaction {0} from library {1}, had a negative rate at 1000 K, 1 bar'.format(rxn,library.label))
            if len(rxn.reactants) == 1 and not rxn.allow_max_rate_violation:
                if k > TST_limit:
                    boo = True
                    logging.error('library reaction {0} from library {1}, exceeds the TST limit at 1000 K, 1 bar of {2} mol/(m3*s) at {3} mol/(m3*s) and did not have allow_max_rate_violation=True'.format(rxn,library.label,TST_limit,k))
            elif len(rxn.reactants) == 2 and not rxn.allow_max_rate_violation:
                if k > collision_limit:
                    boo = True
                    logging.error('library reaction {0} from library {1}, exceeds the collision limit at 1000 K, 1 bar of {2} mol/(m3*s) at {3} mol/(m3*s) and did not have allow_max_rate_violation=True'.format(rxn,library.label,collision_limit,k))
        if boo:
            raise ValueError('library {0} has unreasonable rates'.format(library.label))

    def kinetics_checkReactantAndProductTemplate(self, family_name):
        """
        This test checks whether the reactant and product templates within a family are correctly defined.
        For a reversible family, the reactant and product templates must have matching labels.
        For a non-reversible family, the reactant and product templates must have non-matching labels, otherwise overwriting may occur.
        """
        family = self.database.kinetics.families[family_name]
        if family.ownReverse:
            nose.tools.assert_equal(family.forwardTemplate.reactants, family.forwardTemplate.products)
        else:
            tst = []
            reactant_labels = [reactant.label for reactant in family.forwardTemplate.reactants]
            product_labels = [product.label for product in family.forwardTemplate.products]
            for reactant_label in reactant_labels:
                for product_label in product_labels:
                    tst.append((reactant_label==product_label, "Reactant label {0} matches that of product label {1} in a non-reversible family template.  Please rename product label.".format(reactant_label,product_label)))

            boo = False
            for i in range(len(tst)):
                if tst[i][0]:
                    logging.error(tst[i][1])
                    boo = True

            if boo:
                raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkCdAtomType(self, family_name):
        """
        This test checks that groups containing Cd, CO, CS and Cdd atomtypes are used
        correctly according to their strict definitions
        """
        family = self.database.kinetics.families[family_name]
        targetLabel=['Cd', 'CO', 'CS', 'Cdd']
        targetAtomTypes=[atomTypes[x] for x in targetLabel]

        #ignore product entries that get created from training reactions
        ignore=[]
        if not family.ownReverse:
            for product in family.forwardTemplate.products:
                ignore.append(product)
                ignore.extend(product.children)
        else: ignore=[]
        tst = []
        for entryName, entry in family.groups.entries.items():
            #ignore products
            if entry in ignore: continue
            #ignore LogicOr groups
            if isinstance(entry.item, Group):
                for index, atom in enumerate(entry.item.atoms):
                    for atomtype1 in atom.atomType:
                        if atomtype1 in targetAtomTypes:
                            break
                    #If Cd not found in atomTypes, go to next atom
                    else: continue
                    #Create list of all the atomTypes that should be present in addition or instead of Cd
                    correctAtomList=[]
                    num_of_Dbonds=sum([1 if x.order[0] is 'D' and len(x.order)==1 else 0 for x in atom.bonds.values()])
                    if num_of_Dbonds == 2:
                        correctAtomList.append('Cdd')
                    elif num_of_Dbonds == 1:
                        for ligand, bond in atom.bonds.items():
                            #Ignore ligands that are not double bonded
                            if any([abs(2-order) < 1e-7 for order in bond.order]):
                                for ligAtomType in ligand.atomType:
                                    if ligand.atomType[0].isSpecificCaseOf(atomTypes['O']): correctAtomList.append('CO')
                                    elif ligand.atomType[0].isSpecificCaseOf(atomTypes['S']): correctAtomList.append('CS')

                    #remove duplicates from correctAtom:
                    correctAtomList=list(set(correctAtomList))
                    for correctAtom in correctAtomList:
                        tst.append((atomTypes[correctAtom] in atom.atomType,
                                               """In family {0}, node {1} is missing the atomtype {2} in atom {3} and may be misusing the atomtype Cd, CO, CS, or Cdd.
The following adjList may have atoms in a different ordering than the input file:
{4}
                                            """.format(family_name, entry, correctAtom, index+1, entry.item.toAdjacencyList())))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkUnimolecularGroups(self,family_name):
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

        def getEndFromBackbone(backbone, endLabels):
            """
            :param backbone: :class: Entry for a backbone of molecule
            :param endLabels: Labels in the end groups
            :return: A subgraph representing the end group of the molecule
            """
            #make copy for manipulation
            copyGroup = backbone.item.copy(True)

            #Find the endGroup atoms
            for atom in copyGroup.atoms:
                if atom.label in endLabels:
                    midAtom = atom
                    break

            #find the bonds to break
            bondsToBreak = []
            for atom2, bond in midAtom.bonds.items():
                if atom2.label is None or atom2.label not in endLabels: #
                    bondsToBreak.append(bond)


            for bond in bondsToBreak:
                copyGroup.removeBond(bond)

            #split group into end and backbone fragment
            groups = copyGroup.split()

            #verify group was split correctly and identify the correct end group
            endLabels = set(endLabels)
            for group in groups:
                groupLabels = set(atom.label for atom in group.atoms)
                groupLabels.discard('')
                if endLabels == groupLabels:
                    break
            else:
                raise Exception("Group {0} not split correctly".format(backbone.label))

            return group
        #################################################################################
        family = self.database.kinetics.families[family_name]

        backbone =  family.getBackboneRoots()[0]

        endGroups = family.getEndRoots()

        endLabels = {}
        for endGroup in endGroups:
            labels = []
            for atom in endGroup.item.atoms:
                if atom.label:
                    labels.append(atom.label)
            endLabels[endGroup] = set(labels)

        #get boundary atoms to test that backbones have labels between end groups
        nose.tools.assert_is_not_none(family.boundaryAtoms)

        # set of all end_labels should be backbone label
        backboneLabel = set([])
        for end_label in endLabels.values():
            for label in end_label:
                backboneLabel.add(label)

        #define types of errors
        A = [] #end groups have too many labels
        B = [] #end group lacks necessary label
        C = [] #backbone missing end group labels
        D = [] #backbone missing labels in between groups
        E = [] #backbone tries to define atoms inside end groups
        for group_name, entry in family.groups.entries.items():
            if isinstance(entry.item, Group):
                group = entry.item
                if backbone in family.ancestors(entry):
                    for atom in group.atoms:
                        if atom.label: presentLabels.add(atom.label)
                    #Check C
                    for endGroup, labels in endLabels.items():
                        if not labels.issubset(presentLabels):
                            C.append([endGroup, entry])
                    #check D
                    midAtoms = [group.getLabeledAtom(x)[0] for x in family.boundaryAtoms]
                    pathAtoms = find_shortest_path(midAtoms[0], midAtoms[1])
                    for atom in pathAtoms:
                        if not atom.label:
                            D.append([backbone, entry])
                            break
                    #check E
                    for endGroup, labels in endLabels.items():
                        endFromBackbone = getEndFromBackbone(entry, labels)
                        presentLabels = endFromBackbone.getLabeledAtoms()
                        presentLabels = set(presentLabels.keys())
                        if labels == presentLabels:
                            if not endGroup.item.isIdentical(endFromBackbone):
                                E.append([endGroup, entry])
                        else: raise Exception("Group {0} has split into end group {1}, but does not match any root".format(entry.label, endFromBackbone.toAdjacencyList()))

                else:
                    presentLabels = set([])
                    for endNode, labelledAtoms in endLabels.items():
                        if endNode in family.ancestors(entry):
                            for atom in group.atoms:
                                if atom.label: presentLabels.add(atom.label)
                            #Check A
                            if not presentLabels.issubset(labelledAtoms):
                                A.append([endNode, entry])
                            #Check B
                            if not labelledAtoms.issubset(presentLabels):
                                B.append([endNode, entry])


        #print outputs
        tst = []
        if A != []:
            s = "These end groups have extra labels that their top level end group do not have:"+"\n [root group, error group]"
            for x in A:
                s += '\n'+str(x)
            tst.append((False,s))
        if B != []:
            s = "These end groups are missing labels that their top level end group have:"+"\n [root group, error group]"
            for x in B:
                s += '\n'+str(x)
            tst.append((False,s))
        if C != []:
            s = "These backbone groups are missing labels that are in the end groups:"+"\n [root group, error group]"
            for x in C:
                s += '\n'+str(x)
            tst.append((False,s))
        if D != []:
            s = "These backbone groups are missing labels along the path atoms:"+"\n [root group, error group]"
            for x in D:
                s += '\n'+str(x)
            tst.append((False,s))
        if E != []:
            s = "These backbone have end subgraphs that don't match a root:"+"\n [root group, error group]"
            for x in E:
                s += '\n'+str(x)
            tst.append((False,s))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error occured in databaseTest. Please check log warnings for all error messages.")

    def kinetics_checkSampleDescendsToGroup(self, family_name):
        """
        This test first creates a sample :class:Molecule from a :class:Group. Then it checks
        that this molecule hits the original group or a child when it descends down the tree.
        """
        family = self.database.kinetics.families[family_name]

        tst1 = []
        tst2 = []
        tst3 = []
        #ignore any products
        ignore=[]
        if not family.ownReverse:
            for product in family.forwardTemplate.products:
                ignore.append(product)
                ignore.extend(product.children)
        else: ignore=[]

        #If family is backbone archetype, then we need to merge groups before descending
        roots = family.groups.top
        if len(roots) > len(family.forwardTemplate.reactants):
            backboneRoots = family.getBackboneRoots()
            allBackboneGroups = []
            for backboneRoot in backboneRoots:
                allBackboneGroups.extend(family.getTopLevelGroups(backboneRoot))
            #list of numbered of labelled atoms for allBackboneGroups
            backboneSizes = [len(backbone.item.getLabeledAtoms()) for backbone in allBackboneGroups]

            #pick a backbone that is two labelled atoms larger than the smallest
            if min(backboneSizes) + 2 in backboneSizes:
                backboneSample = allBackboneGroups[backboneSizes.index(min(backboneSizes) + 2)]
            #or if it doesn't exist, pick the largest backbone
            else:
                backboneSample = allBackboneGroups[backboneSizes.index(max(backboneSizes))]
            mergesNecessary = True
        else: mergesNecessary = False

        #If atom has too many benzene rings, we currently have trouble making sample atoms
        skipped = []
        for entryName, entry in family.groups.entries.items():
            if entry in ignore: continue
            elif isinstance(entry.item, Group):
                ancestors=family.ancestors(entry)
                if ancestors: root = ancestors[-1] #top level root will be last one in ancestors
                else: root = entry
                try:
                    if mergesNecessary and root not in backboneRoots: #we may need to merge
                        mergedGroup = backboneSample.item.mergeGroups(entry.item)
                        sampleMolecule = mergedGroup.makeSampleMolecule()
                    else:
                        sampleMolecule = entry.item.makeSampleMolecule()

                    #test accessibility here
                    atoms = sampleMolecule.getLabeledAtoms()
                    match = family.groups.descendTree(sampleMolecule, atoms, strict=True, root = root)
                    tst1.append((match, "Group {0} does not match its root node, {1}".format(entryName, root.label)))
                    if tst1[-1][0] is not None:
                        tst2.append((entry, [match]+family.groups.ancestors(match), """In group {0}, a sample molecule made from node {1} returns node {2} when descending the tree.
Sample molecule AdjList:
{3}

Origin Group AdjList:
{4}{5}{6}

Matched group AdjList:
{7}
        """.format(family_name,
                   entry.label,
                   match.label,
                   sampleMolecule.toAdjacencyList(),
                   entry.item.toAdjacencyList(),
                   "\n\nBackbone Group Adjlist:\n" + backboneSample.label +'\n' if mergesNecessary and root not in backboneRoots else '',
                   backboneSample.item.toAdjacencyList() if mergesNecessary and root not in backboneRoots else '',
                   match.item.toAdjacencyList())))

                except UnexpectedChargeError as e:
                    tst3.append((False, """In family {0}, a sample molecule made from node {1} returns an unexpectedly charged molecule:
Sample molecule AdjList:
{2}

Origin Group AdjList:
{3}{4}{5}""".format(family_name,
                    entry.label,
                    e.graph.toAdjacencyList(),
                    entry.item.toAdjacencyList(),
                    "\n\nBackbone Group Adjlist:\n" + backboneSample.label +'\n' if mergesNecessary and root not in backboneRoots else '',
                    backboneSample.item.toAdjacencyList() if mergesNecessary and root not in backboneRoots else '')
                   ))

                except ImplicitBenzeneError:
                    skipped.append(entryName)

        #print out entries skipped from exception we can't currently handle
        if skipped:
            print("These entries were skipped because too big benzene rings or has nitrogen sample atom:")
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

    def general_checkNodesFoundInTree(self, group_name, group):
        """
        This test checks whether nodes are found in the tree, with proper parents.
        """
        for nodeName, nodeGroup in group.entries.items():
            ascendParent = nodeGroup
            # Check whether the node has proper parents unless it is the top reactant or product node
            tst1 = []
            tst2 = []
            tst3 = []
            while ascendParent not in group.top:
                child = ascendParent
                ascendParent = ascendParent.parent
                tst1.append((ascendParent is not None, "Node {node} in {group} group was found in the tree without a proper parent.".format(node=child, group=group_name)))
                if tst1[-1] is not None:
                    tst2.append((child in ascendParent.children, "Node {node} in {group} group was found in the tree without a proper parent.".format(node=nodeName, group=group_name)))
                    tst3.append((child is ascendParent, "Node {node} in {group} is a parent to itself".format(node=nodeName, group=group_name)))

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

    def general_checkGroupsNonidentical(self, group_name, group):
        """
        This test checks whether nodes found in the group are nonidentical.
        """
        entriesCopy = copy(group.entries)
        tst = []
        for nodeName, nodeGroup in group.entries.items():
            del entriesCopy[nodeName]
            for nodeNameOther, nodeGroupOther in entriesCopy.items():
                group.matchNodeToNode(nodeGroup,nodeGroupOther)
                tst.append((group.matchNodeToNode(nodeGroup, nodeGroupOther), "Node {node} in {group} group was found to be identical to node {nodeOther}".format(node=nodeName, group=group_name, nodeOther=nodeNameOther)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_checkChildParentRelationships(self, group_name, group):
        """
        This test checks that nodes' parent-child relationships are correct in the database.
        """
        tst1 = []
        tst2 = []
        for nodeName, childNode in group.entries.items():
            #top nodes and product nodes don't have parents by definition, so they get an automatic pass:
            if childNode in group.top: continue
            parentNode = childNode.parent
            # Check whether the node has proper parents unless it is the top reactant or product node
            # The parent should be more general than the child
            tst1.append((group.matchNodeToChild(parentNode, childNode),
                            "In {group} group, node {parent} is not a proper parent of its child {child}.".format(group=group_name, parent=parentNode, child=nodeName))
)

            #check that parentNodes which are LogicOr do not have an ancestor that is a Group
            #If it does, then the childNode must also be a child of the ancestor
            if isinstance(parentNode.item, LogicOr):
                ancestorNode = parentNode
                while ancestorNode not in group.top and isinstance(ancestorNode.item, LogicOr):
                    ancestorNode = ancestorNode.parent
                if isinstance(ancestorNode.item, Group) and tst1[-1][0]:
                    tst2.append((group.matchNodeToChild(ancestorNode, childNode),
                                    "In {group} group, node {ancestor} is not a proper ancestor of its child {child}.".format(group=group_name, ancestor=ancestorNode, child=nodeName))
)

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

    def general_checkSiblingsForParents(self, group_name, group):
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
                for child2 in node.children[index+1:]:
                    tst.append((group.matchNodeToChild(child1, child2),
                                            "In {0} group, node {1} is a parent of {2}, but they are written as siblings.".format(group_name, child1, child2)))

        boo = False
        for i in range(len(tst)):
            if tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_checkCdAtomType(self, group_name, group):
        """
        This test checks that groups containing Cd, CO, CS and Cdd atomtypes are used
        correctly according to their strict definitions
        """
        targetLabel=['Cd', 'CO', 'CS', 'Cdd']
        targetAtomTypes=[atomTypes[x] for x in targetLabel]
        tst = []
        for entryName, entry in group.entries.items():
            if isinstance(entry.item, Group):
                for index, atom in enumerate(entry.item.atoms):
                    for atomtype1 in atom.atomType:
                        if atomtype1 in targetAtomTypes:
                            break
                    #If Cd not found in atomTypes, go to next atom
                    else: continue
                    #figure out what the correct atomType is
                    correctAtomList=[]
                    num_of_Dbonds=sum([1 if x.order[0] is 'D' and len(x.order)==1 else 0 for x in atom.bonds.values()])
                    if num_of_Dbonds == 2:
                        correctAtomList.append('Cdd')
                    elif num_of_Dbonds == 1:
                        for ligand, bond in atom.bonds.items():
                            #Ignore ligands that are not double bonded
                            if any([abs(2-order) < 1e-7 for order in bond.order]):
                                for ligAtomType in ligand.atomType:
                                    if ligand.atomType[0].isSpecificCaseOf(atomTypes['O']): correctAtomList.append('CO')
                                    elif ligand.atomType[0].isSpecificCaseOf(atomTypes['S']): correctAtomList.append('CS')

                    #remove duplicates from correctAtom:
                    correctAtomList=list(set(correctAtomList))
                    for correctAtom in correctAtomList:
                        tst.append((atomTypes[correctAtom] in atom.atomType,
                                                """In group {0}, node {1} is missing the atomtype {2} in atom {3} and may be misusing the atomtype Cd, CO, CS, or Cdd.
The following adjList may have atoms in a different ordering than the input file:
{4}
                                            """.format(group_name, entry, correctAtom, index+1, entry.item.toAdjacencyList())))

        boo = False
        for i in range(len(tst)):
            if not tst[i][0]:
                logging.error(tst[i][1])
                boo = True

        if boo:
            raise ValueError("Error Occurred")

    def general_checkSampleDescendsToGroup(self, group_name, group):
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
                        sampleMolecule = entry.item.makeSampleMolecule()
                    except:
                        logging.error("Problem making sample molecule for group {}\n{}".format(entryName, entry.item.toAdjacencyList()))
                        raise
                    #for now ignore sample atoms that use nitrogen types
                    nitrogen = False
                    for atom in sampleMolecule.atoms:
                        if atom.isNitrogen(): nitrogen = True
                    if nitrogen:
                        skipped.append(entryName)
                        continue

                    atoms = sampleMolecule.getLabeledAtoms()
                    match = group.descendTree(sampleMolecule, atoms, strict=True)
                    tst1.append((match, "Group {0} does not match its root node, {1}".format(entryName, group.top[0])))
                    tst2.append((entry, [match]+group.ancestors(match), """In group {0}, a sample molecule made from node {1} returns node {2} when descending the tree.
Sample molecule AdjList:
{3}

Origin Group AdjList:
{4}

Matched group AdjList:
{5}
""".format(group_name,
           entry,
           match,
           sampleMolecule.toAdjacencyList(),
           entry.item.toAdjacencyList(),
           match.item.toAdjacencyList())))

            except UnexpectedChargeError as e:
                tst3.append((False, """In family {0}, a sample molecule made from node {1} returns an unexpectedly charged molecule:
Sample molecule AdjList:
{2}

Origin Group AdjList:
{3}""".format(group_name,
                    entry.label,
                    e.graph.toAdjacencyList(),
                    entry.item.toAdjacencyList())))

            except ImplicitBenzeneError:
                skipped.append(entryName)

        #print out entries skipped from exception we can't currently handle
        if skipped:
            print("These entries were skipped because too big benzene rings or has nitrogen sample atom:")
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
