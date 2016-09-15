"""
This scripts runs tests on the database 
"""
import os.path
import logging
from external.wip import work_in_progress
from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from copy import copy, deepcopy
from rmgpy.data.base import LogicOr
from rmgpy.molecule import Group
from rmgpy.molecule.atomtype import atomTypes

import nose
import nose.tools


class TestDatabase():  # cannot inherit from unittest.TestCase if we want to use nose test generators
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
        for family_name, family in self.database.kinetics.families.iteritems():

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
            
            for depository in family.depositories:
                
                test = lambda x: self.kinetics_checkAdjlistsNonidentical(depository)
                test_name = "Kinetics {1} Depository: check adjacency lists are nonidentical?".format(family_name, depository.label)
                test.description = test_name
                self.compat_func_name = test_name
                yield test, depository.label
        
        for library_name, library in self.database.kinetics.libraries.iteritems():
            
            test = lambda x: self.kinetics_checkAdjlistsNonidentical(library)
            test_name = "Kinetics library {0}: check adjacency lists are nonidentical?".format(library_name)
            test.description = test_name
            self.compat_func_name = test_name
            yield test, library_name
        
    def test_thermo(self):
        for group_name, group in self.database.thermo.groups.iteritems():
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

    def test_solvation(self):
        for group_name, group in self.database.solvation.groups.iteritems():
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

    def test_statmech(self):
        for group_name, group in self.database.statmech.groups.iteritems():
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

    def test_transport(self):
        for group_name, group in self.database.transport.groups.iteritems():
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
            
    # These are the actual tests, that don't start with a "test_" name:
    def kinetics_checkCorrectNumberofNodesInRules(self, family_name):
        """
        This test ensures that each rate rule contains the proper number of nodes according to the family it originates.
        """
        family = self.database.kinetics.families[family_name]
        expectedNumberNodes = len(family.getRootTemplate())
        for label, entries in family.rules.entries.iteritems():
            for entry in entries:
                nodes = label.split(';')
                nose.tools.assert_equal(len(nodes), expectedNumberNodes, "Wrong number of groups or semicolons in family {family} rule {entry}.  Should be {num_nodes}".format(family=family_name, entry=entry, num_nodes=expectedNumberNodes))

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
        
        for label, entries in family.rules.entries.iteritems():
            for entry in entries:
                nodes = label.split(';')
                for i, node in enumerate(nodes):
                    nose.tools.assert_true(node in family.groups.entries, "In {family} family, no group definition found for label {label} in rule {entry}".format(family=family_name, label=node, entry=entry))
                    nose.tools.assert_true(family.groups.entries[node] in topDescendants[i], "In {family} family, rule {entry} was found with groups out of order.  The correct order for a rule should be subgroups of {top}.".format(family=family_name, entry=entry, top=topGroupOrder))
                                        
    def kinetics_checkGroupsFoundInTree(self, family_name):
        """
        This test checks whether groups are found in the tree, with proper parents.
        """
        family = self.database.kinetics.families[family_name]
        for nodeName, nodeGroup in family.groups.entries.iteritems():
            nose.tools.assert_false('[' in nodeName or ']' in nodeName, "Group {group} in {family} family contains square brackets [ ] in the label, which are not allowed.".format(group=nodeName, family=family_name))
            ascendParent = nodeGroup
            # Check whether the node has proper parents unless it is the top reactant or product node
            while ascendParent not in family.groups.top and ascendParent not in family.forwardTemplate.products:
                child = ascendParent
                ascendParent = ascendParent.parent
                nose.tools.assert_true(ascendParent is not None, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=child, family=family_name))
                nose.tools.assert_true(child in ascendParent.children, "Group {group} in {family} family was found in the tree without a proper parent.".format(group=nodeName, family=family_name))
                nose.tools.assert_false(child is ascendParent, "Group {group} in {family} family is a parent to itself".format(group=nodeName, family=family_name))

    def kinetics_checkGroupsNonidentical(self, family_name):
        """
        This test checks that the groups are non-identical.
        """
        from rmgpy.data.base import Database
        originalFamily = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = originalFamily.groups.entries
        entriesCopy = copy(family.entries)
        for nodeName, nodeGroup in family.entries.iteritems():
            del entriesCopy[nodeName]
            for nodeNameOther, nodeGroupOther in entriesCopy.iteritems():
                nose.tools.assert_false(family.matchNodeToNode(nodeGroup, nodeGroupOther), "Group {group} in {family} family was found to be identical to group {groupOther}".format(group=nodeName, family=family_name, groupOther=nodeNameOther))

    def kinetics_checkChildParentRelationships(self, family_name):
        """
        This test checks that groups' parent-child relationships are correct in the database.
        """
        from rmgpy.data.base import Database
        originalFamily = self.database.kinetics.families[family_name]
        family = Database()
        family.entries = originalFamily.groups.entries
        for nodeName, childNode in family.entries.iteritems():
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
            nose.tools.assert_true(family.matchNodeToChild(parentNode, childNode),
                            "In {family} family, group {parent} is not a proper parent of its child {child}.".format(family=family_name, parent=parentNode, child=nodeName))

            #check that parentNodes which are LogicOr do not have an ancestor that is a Group
            #If it does, then the childNode must also be a child of the ancestor
            if isinstance(parentNode.item, LogicOr):
                ancestorNode = parentNode
                while ancestorNode not in originalFamily.groups.top and isinstance(ancestorNode.item, LogicOr):
                    ancestorNode = ancestorNode.parent
                if isinstance(ancestorNode.item, Group):
                    nose.tools.assert_true(family.matchNodeToChild(ancestorNode, childNode),
                                    "In {family} family, group {ancestor} is not a proper ancestor of its child {child}.".format(family=family_name, ancestor=ancestorNode, child=nodeName))

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
        for nodeName, node in family.entries.iteritems():
            #Some families also construct a 2-level trees for the products
            #(root with all entries down one level) We don't care about this
            #tree as it is not used in searching, so we ignore products
            if node in originalFamily.forwardTemplate.products: continue
            for index, child1 in enumerate(node.children):
                for child2 in node.children[index+1:]:
                    nose.tools.assert_false(family.matchNodeToChild(child1, child2),
                                            "In family {0}, node {1} is a parent of {2}, but they are written as siblings.".format(family_name, child1, child2))

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
                    
        # Go through all species to make sure they are nonidentical
        speciesList = speciesDict.values()
        labeledAtoms = [species.molecule[0].getLabeledAtoms() for species in speciesList]
        for i in range(len(speciesList)):
            for j in range(i+1,len(speciesList)):
                    initialMap = {}
                    try:
                        for atomLabel in labeledAtoms[i]:
                            initialMap[labeledAtoms[i][atomLabel]] = labeledAtoms[j][atomLabel]
                    except KeyError:
                        # atom labels did not match, therefore not a match
                        continue
                    
                    nose.tools.assert_false(speciesList[i].molecule[0].isIsomorphic(speciesList[j].molecule[0], initialMap), "Species {0} and species {1} in {2} database were found to be identical.".format(speciesList[i].label,speciesList[j].label,database.label))

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
            reactant_labels = [reactant.label for reactant in family.forwardTemplate.reactants]
            product_labels = [product.label for product in family.forwardTemplate.products]
            for reactant_label in reactant_labels:
                for product_label in product_labels:
                    nose.tools.assert_false(reactant_label==product_label, "Reactant label {0} matches that of product label {1} in a non-reversible family template.  Please rename product label.".format(reactant_label,product_label))
        
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

        for entryName, entry in family.groups.entries.iteritems():
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
                        for ligand, bond in atom.bonds.iteritems():
                            #Ignore ligands that are not double bonded
                            if 'D' in bond.order:
                                for ligAtomType in ligand.atomType:
                                    if ligand.atomType[0].isSpecificCaseOf(atomTypes['O']): correctAtomList.append('CO')
                                    elif ligand.atomType[0].isSpecificCaseOf(atomTypes['S']): correctAtomList.append('CS')

                    #remove duplicates from correctAtom:
                    correctAtomList=list(set(correctAtomList))
                    for correctAtom in correctAtomList:
                        nose.tools.assert_true(atomTypes[correctAtom] in atom.atomType,
                                               """In family {0}, node {1} is missing the atomtype {2} in atom {3} and may be misusing the atomtype Cd, CO, CS, or Cdd.
The following adjList may have atoms in a different ordering than the input file:
{4}
                                            """.format(family_name, entry, correctAtom, index+1, entry.item.toAdjacencyList()))


    def general_checkNodesFoundInTree(self, group_name, group):
        """
        This test checks whether nodes are found in the tree, with proper parents.
        """
        for nodeName, nodeGroup in group.entries.iteritems():
            ascendParent = nodeGroup
            # Check whether the node has proper parents unless it is the top reactant or product node
            while ascendParent not in group.top:
                child = ascendParent
                ascendParent = ascendParent.parent
                nose.tools.assert_true(ascendParent is not None, "Node {node} in {group} group was found in the tree without a proper parent.".format(node=child, group=group_name))
                nose.tools.assert_true(child in ascendParent.children, "Node {node} in {group} group was found in the tree without a proper parent.".format(node=nodeName, group=group_name))
                nose.tools.assert_false(child is ascendParent, "Node {node} in {group} is a parent to itself".format(node=nodeName, group=group_name))
    
    def general_checkGroupsNonidentical(self, group_name, group):
        """
        This test checks whether nodes found in the group are nonidentical.
        """
        entriesCopy = copy(group.entries)
        for nodeName, nodeGroup in group.entries.iteritems():
            del entriesCopy[nodeName]
            for nodeNameOther, nodeGroupOther in entriesCopy.iteritems():
                try: 
                    group.matchNodeToNode(nodeGroup,nodeGroupOther)
                except:
                    print nodeName
                    print nodeNameOther
                    pass
                nose.tools.assert_false(group.matchNodeToNode(nodeGroup, nodeGroupOther), "Node {node} in {group} group was found to be identical to node {nodeOther}".format(node=nodeName, group=group_name, nodeOther=nodeNameOther))
    
    def general_checkChildParentRelationships(self, group_name, group):
        """
        This test checks that nodes' parent-child relationships are correct in the database.
        """
        for nodeName, childNode in group.entries.iteritems():
            #top nodes and product nodes don't have parents by definition, so they get an automatic pass:
            if childNode in group.top: continue
            parentNode = childNode.parent
            # Check whether the node has proper parents unless it is the top reactant or product node
            # The parent should be more general than the child
            nose.tools.assert_true(group.matchNodeToChild(parentNode, childNode),
                            "In {group} group, node {parent} is not a proper parent of its child {child}.".format(group=group_name, parent=parentNode, child=nodeName))

            #check that parentNodes which are LogicOr do not have an ancestor that is a Group
            #If it does, then the childNode must also be a child of the ancestor
            if isinstance(parentNode.item, LogicOr):
                ancestorNode = parentNode
                while ancestorNode not in group.top and isinstance(ancestorNode.item, LogicOr):
                    ancestorNode = ancestorNode.parent
                if isinstance(ancestorNode.item, Group):
                    nose.tools.assert_true(group.matchNodeToChild(ancestorNode, childNode),
                                    "In {group} group, node {ancestor} is not a proper ancestor of its child {child}.".format(group=group_name, ancestor=ancestorNode, child=nodeName))

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
        for nodeName, node in group.entries.iteritems():
            for index, child1 in enumerate(node.children):
                for child2 in node.children[index+1:]:
                    nose.tools.assert_false(group.matchNodeToChild(child1, child2),
                                            "In {0} group, node {1} is a parent of {2}, but they are written as siblings.".format(group_name, child1, child2))

    def general_checkCdAtomType(self, group_name, group):
        """
        This test checks that groups containing Cd, CO, CS and Cdd atomtypes are used
        correctly according to their strict definitions
        """
        targetLabel=['Cd', 'CO', 'CS', 'Cdd']
        targetAtomTypes=[atomTypes[x] for x in targetLabel]

        for entryName, entry in group.entries.iteritems():
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
                        for ligand, bond in atom.bonds.iteritems():
                            #Ignore ligands that are not double bonded
                            if 'D' in bond.order:
                                for ligAtomType in ligand.atomType:
                                    if ligand.atomType[0].isSpecificCaseOf(atomTypes['O']): correctAtomList.append('CO')
                                    elif ligand.atomType[0].isSpecificCaseOf(atomTypes['S']): correctAtomList.append('CS')

                    #remove duplicates from correctAtom:
                    correctAtomList=list(set(correctAtomList))
                    for correctAtom in correctAtomList:
                        nose.tools.assert_true(atomTypes[correctAtom] in atom.atomType,
                                                """In group {0}, node {1} is missing the atomtype {2} in atom {3} and may be misusing the atomtype Cd, CO, CS, or Cdd.
The following adjList may have atoms in a different ordering than the input file:
{4}
                                            """.format(group_name, entry, correctAtom, index+1, entry.item.toAdjacencyList()))

if __name__ == '__main__':
    nose.run(argv=[__file__, '-v', '--nologcapture'], defaultTest=__name__)
