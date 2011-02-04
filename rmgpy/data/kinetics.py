#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2010 Prof. William H. Green (whgreen@mit.edu) and the
#   RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

import os
import os.path
import logging

from base import *

from chempy.reaction import Reaction, ReactionError
from chempy.kinetics import *
from chempy.pattern import BondPattern, MoleculePattern
from chempy.molecule import Bond

################################################################################

class UndeterminableKineticsError(ReactionError):
    """
    An exception raised when attempts to estimate appropriate kinetic parameters
    for a chemical reaction are unsuccessful.
    """
    def __init__(self, reaction, message=''):
        new_message = 'Kinetics could not be determined. '+message
        ReactionError.__init__(self,reaction,new_message)

class InvalidActionError(Exception):
    """
    An exception to be raised when an invalid action is encountered in a
    reaction recipe.
    """
    pass

################################################################################

class ReactionRecipe:
    """
    Represent a list of actions that, when executed, result in the conversion
    of a set of reactants to a set of products. There are currently five such
    actions:

    ============= ============================= ================================
    Action Name   Arguments                     Description
    ============= ============================= ================================
    CHANGE_BOND   `center1`, `order`, `center2` change the bond order of the bond between `center1` and `center2` by `order`; do not break or form bonds
    FORM_BOND     `center1`, `order`, `center2` form a new bond between `center1` and `center2` of type `order`
    BREAK_BOND    `center1`, `order`, `center2` break the bond between `center1` and `center2`, which should be of type `order`
    GAIN_RADICAL  `center`, `radical`           increase the number of free electrons on `center` by `radical`
    LOSE_RADICAL  `center`, `radical`           decrease the number of free electrons on `center` by `radical`
    ============= ============================= ================================

    The actions are stored as a list in the `actions` attribute. Each action is
    a list of items; the first is the action name, while the rest are the
    action parameters as indicated above.
    """

    def __init__(self, actions=None):
        self.actions = actions or []

    def addAction(self, action):
        """
        Add an `action` to the reaction recipe, where `action` is a list
        containing the action name and the required parameters, as indicated in
        the table above.
        """
        self.actions.append(action)

    def getReverse(self):
        """
        Generate a reaction recipe that, when applied, does the opposite of
        what the current recipe does, i.e., it is the recipe for the reverse
        of the reaction that this is the recipe for.
        """
        other = ReactionRecipe()
        for action in self.actions:
            if action[0] == 'CHANGE_BOND':
                other.addAction(['CHANGE_BOND', action[1], str(-int(action[2])), action[3]])
            elif action[0] == 'FORM_BOND':
                other.addAction(['BREAK_BOND', action[1], action[2], action[3]])
            elif action[0] == 'BREAK_BOND':
                other.addAction(['FORM_BOND', action[1], action[2], action[3]])
            elif action[0] == 'LOSE_RADICAL':
                other.addAction(['GAIN_RADICAL', action[1], action[2]])
            elif action[0] == 'GAIN_RADICAL':
                other.addAction(['LOSE_RADICAL', action[1], action[2]])
        return other

    def __apply(self, struct, doForward, unique):
        """
        Apply the reaction recipe to the set of molecules contained in
        `structure`, a single Structure object that contains one or more
        structures. The `doForward` parameter is used to indicate
        whether the forward or reverse recipe should be applied. The atoms in
        the structure should be labeled with the appropriate atom centers.
        """

        pattern = isinstance(struct, MoleculePattern)

        for action in self.actions:
            if action[0] in ['CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND']:

                # We are about to change the connectivity of the atoms in
                # struct, which invalidates any existing vertex connectivity
                # information; thus we reset it
                struct.resetConnectivityValues()

                label1, info, label2 = action[1:]
                
                # Find associated atoms
                atom1 = struct.getLabeledAtom(label1)
                atom2 = struct.getLabeledAtom(label2)
                if atom1 is None or atom2 is None or atom1 is atom2:
                    raise InvalidActionError('Invalid atom labels encountered.')
                
                # Apply the action
                if action[0] == 'CHANGE_BOND':
                    info = int(info)
                    bond = struct.getBond(atom1, atom2)
                    if doForward:
                        atom1.applyAction(['CHANGE_BOND', label1, info, label2])
                        atom2.applyAction(['CHANGE_BOND', label1, info, label2])
                        bond.applyAction(['CHANGE_BOND', label1, info, label2])
                    else:
                        atom1.applyAction(['CHANGE_BOND', label1, -info, label2])
                        atom2.applyAction(['CHANGE_BOND', label1, -info, label2])
                        bond.applyAction(['CHANGE_BOND', label1, -info, label2])
                elif (action[0] == 'FORM_BOND' and doForward) or (action[0] == 'BREAK_BOND' and not doForward):
                    bond = BondPattern(order=['S']) if pattern else Bond(order='S')
                    struct.addBond(atom1, atom2, bond)
                    atom1.applyAction(['FORM_BOND', label1, info, label2])
                    atom2.applyAction(['FORM_BOND', label1, info, label2])
                elif (action[0] == 'BREAK_BOND' and doForward) or (action[0] == 'FORM_BOND' and not doForward):
                    if not struct.hasBond(atom1, atom2):
                        raise InvalidActionError('Attempted to remove a nonexistent bond.')
                    bond = struct.getBond(atom1, atom2)
                    struct.removeBond(atom1, atom2)
                    atom1.applyAction(['BREAK_BOND', label1, info, label2])
                    atom2.applyAction(['BREAK_BOND', label1, info, label2])

            elif action[0] in ['LOSE_RADICAL', 'GAIN_RADICAL']:

                label, change = action[1:]
                change = int(change)
                
                # Find associated atom
                atom = struct.getLabeledAtom(label)
                if atom is None:
                    raise InvalidActionError('Unable to find atom with label "%s" while applying reaction recipe.' % label)
                
                # Apply the action
                for i in range(change):
                    if (action[0] == 'GAIN_RADICAL' and doForward) or (action[0] == 'LOSE_RADICAL' and not doForward):
                        atom.applyAction(['GAIN_RADICAL', label, 1])
                    elif (action[0] == 'LOSE_RADICAL' and doForward) or (action[0] == 'GAIN_RADICAL' and not doForward):
                        atom.applyAction(['LOSE_RADICAL', label, 1])

            else:
                raise InvalidActionError('Unknown action "' + action[0] + '" encountered.')

        struct.updateConnectivityValues()

    def applyForward(self, struct, unique=True):
        """
        Apply the forward reaction recipe to `molecule`, a single
        :class:`Molecule` object.
        """
        return self.__apply(struct, True, unique)

    def applyReverse(self, struct, unique=True):
        """
        Apply the reverse reaction recipe to `molecule`, a single
        :class:`Molecule` object.
        """
        return self.__apply(struct, False, unique)

################################################################################

class KineticsEntry(DataEntry):
    """
    A single entry in the kinetics database. Each entry either contains
    a kinetics `model` or a string label of a `node` to look at for
    kinetics information.
    """

    def __init__(self, model=None, node='', index=0, label='', shortComment='', longComment='', history=None):
        DataEntry.__init__(self, index, label, shortComment, longComment, history)
        self.model = model
        self.node = node

################################################################################

class ReactionFamily(Database):
    """
    Represent a reaction family: a set of reactions with similar chemistry, and
    therefore similar reaction rates. Besides the dictionary, tree, and library
    inherited from :class:`Database`, the attributes are:

    =================== ======================= ================================
    Attribute           Type                    Description
    =================== ======================= ================================
    `label`             ``str``                 The name of the reaction family
    `forwardTemplate`   :class:`Reaction`       The forward reaction template
    `forwardRecipe`     :class:`ReactionRecipe` The steps to take when applying the forward reaction to a set of reactants
    `reverseTemplate`   :class:`Reaction`       The reverse reaction template
    `reverseRecipe`     :class:`ReactionRecipe` The steps to take when applying the reverse reaction to a set of reactants
    `forbidden`         ``dict``                (Optional) Forbidden product structures in either direction
    =================== ======================= ================================

    There are a few reaction families that are their own reverse (hydrogen
    abstraction and intramolecular hydrogen migration); for these
    `reverseTemplate` and `reverseRecipe` will both be ``None``.
    """

    def __init__(self, label=''):
        Database.__init__(self)
        # calling  Database.__init__(self) sets:
        # 	self.dictionary = Dictionary()
        # 	self.library = Library()
        # 	self.tree = Tree()
        self.label = label
        self.forwardTemplate = None
        self.forwardRecipe = None
        self.reverseTemplate = None
        self.reverseRecipe = None
        self.forbidden = None

    def __str__(self):
        return '<ReactionFamily "%s">' % (self.label)

    def load(self, path):
        """
        Load a reaction family located in the directory `path`. The family
        consists of the files ``dictionary.txt``, ``tree.txt``, 
        ``rateLibrary.txt``, ``reactionAdjList.txt``, and optionally
        ``forbiddenGroups.txt``.
        """
        # Generate paths to files in the database
        dictstr = os.path.join(path, 'dictionary.txt')
        treestr = os.path.join(path, 'tree.txt')
        libstr  = os.path.join(path, 'rateLibrary.txt')
        tempstr = os.path.join(path, 'reactionAdjList.txt')
        forbstr = os.path.join(path, 'forbiddenGroups.txt')

        #: The path of the database that was loaded.
        self._path = path

        # Load the dictionary and tree using the generic methods
        # We can't use the generic method to load the library because it has
        # the type ('Arrhenius_EP') as the first meaningful line
        Database.load(self, dictstr, treestr, '')

        # Load the forbidden groups if the file 'forbiddenGroups.txt' is present
        # This file has the form of a standard dictionary so we can use the
        # standard dictionary loading function
        if os.path.exists(forbstr):
            self.forbidden = Dictionary()
            self.forbidden.load(forbstr)
        
        # Load the reaction template information and generate the reverse family
        # This requires that the dictionary and tree be loaded
        self.loadTemplate(tempstr)

        # Process the data in the library
        lines = self.library.load(libstr)

        # pop off the first line and check that it's 'Arrhenius_EP'
        line = lines.pop(0)
        if line != 'Arrhenius_EP':
            raise InvalidDatabaseError("Was expecting 'Arrhenius_EP' as first line, but got %s, in %s"%(line,libstr))

        #figure out how many labels there are
        test_line = lines[0].split()
        for token_no, token in enumerate(test_line):
            # skip token_no=0 because it's the label (but may match the regular expression)
            if token_no and re.match('^[0-9\-.]*$',token):
                # found the Temperature range at token_no
                number_of_groups = token_no-1
                logging.log(0, "Deduced there are %d groups %s in %s"%(number_of_groups,test_line[1:token_no],libstr))
                break
        else: # didn't break
            raise InvalidDatabaseError("Unable to figure out how many groups in %s using line %s"%(libstr,' '.join(test_line)))

        self.library.parse(lines, number_of_groups)
        self.processLibraryData()

        # Check for well-formedness
        if not self.isWellFormed():
            raise InvalidDatabaseError('Database at "%s" is not well-formed.' % (path))

    def processLibraryData(self):
        """
        Convert the data in the library from a string/unicode object to either
        an :class:`ArrheniusEPModel` object or a list of [link, comment]
        string pairs. This function is generally called in the course of
        loading a database from files.
        """

        numReactants = len(self.forwardTemplate.reactants)

        for label, item in self.library.iteritems():

            if item is None:
                pass
            elif not item.__class__ is tuple:
                raise InvalidDatabaseError('Kinetics library should be tuple at this point. Instead got %r'%data)
            else:
                index,data = item # break apart tuple, recover the 'index' - the beginning of the line in the library file.
                # Is't it dangerous having a local variable with the same name as a module?
                # what if we want to raise another InvalidDatabaseError() ?
                if not ( data.__class__ is str or data.__class__ is unicode) :
                    raise InvalidDatabaseError('Kinetics library data format is unrecognized.')

                items = data.split()
                try:
                    kineticData = [];
                    # First item is temperature range
                    kineticData.extend(items[0].split('-'))
                    if len(kineticData) == 2:
                        kineticData[0] = float(kineticData[0])
                        kineticData[1] = float(kineticData[1])
                    elif len(kineticData) == 1:
                        kineticData = [float(items[0]), float(items[0])]
                    # Middle items are Arrhenius + Evans-Polanyi data
                    for i in range(1, 5):
                        kineticData.append(float(items[i]))
                    for i in range(5, 9):
                        # Convert multiplicative uncertainties to additive
                        # uncertainties if needed
                        if items[i][0] == '*':
                            kineticData.append((float(items[i][1:]) - 1.0) * float(items[i-4]))
                        else:
                            kineticData.append(float(items[i]))
                    # Final item before comment is quality
                    kineticData.append(int(items[9]))
                    # Everything else is a comment
                    comment = ' '.join(items[10:]).replace('"', '').strip()

                    # Convert data to ArrheniusEPModel object
                    if len(kineticData) != 11:
                        raise Exception('Invalid list of kinetic data. Should be a list of numbers of length 11; instead got %s'%data)
                    Tmin, Tmax, A, n, alpha, E0, dA, dn, dalpha, dE0, rank = kineticData
                    # Get units of preexponential
                    originalUnits = 's^-1'; desiredUnits = 's^-1'
                    if numReactants == 2:
                        originalUnits = 'cm^3/(mol*s)'
                    elif numReactants > 2:
                        originalUnits = 'cm^%s/(mol^%s*s)' % ((numReactants-1)*3, numReactants-1)
                    # Convert parameters to proper units
                    Tmin = float(pq.Quantity(Tmin, 'K').simplified)
                    Tmax = float(pq.Quantity(Tmax, 'K').simplified)
                    A = float(pq.Quantity(A, originalUnits).simplified)
                    E0 = float(pq.Quantity(E0, 'kcal/mol').simplified)
                    n = float(pq.Quantity(n, '').simplified)
                    alpha = float(pq.Quantity(alpha, '').simplified)
                    # Construct ArrheniusEPModel object
                    kinetics = ArrheniusEPModel(A=A, n=n, alpha=alpha, E0=E0)
                    kinetics.Tmin = Tmin
                    kinetics.Tmax = Tmax
                    kinetics.comment = comment
                    
                    self.library[label] = KineticsEntry(
                        index=int(index),
                        label=label,
                        model=kinetics,
                        shortComment = comment
                    )
                    
                except (ValueError, IndexError), e:
                    # Data represents a link to a different node that contains
                    # the data to use
                    database.library[label] = KineticsEntry(
                        node=items[0],
                        index=int(index),
                        label=label,
                        shortComment=item[len(items[0])+1:].replace('"', '').strip(),
                    )

    def loadTemplate(self, path):
        """
        Load and process a reaction template file located at `path`. This file
        is part of every reaction family.
        """

        # Process the template file, removing comments and empty lines
        info = ''
        frec = None
        try:
            frec = open(path, 'r')
            for line in frec:
                line = removeCommentFromLine(line).strip()
                if len(line) > 0:
                    info += line + '\n'
        except InvalidDatabaseError, e:
            logging.exception(str(e))
            return
        except IOError, e:
            logging.exception('Database file "' + e.filename + '" not found.')
            raise
        finally:
            if frec: frec.close()

        lines = info.splitlines()

        # First line is reaction template, of the form A + B -> C + D
        reactants = []; products = []; species = []
        atArrow = False
        items = lines[0].split(); items.extend('+')
        for item in items:
            if item == '+' or item == '->':
                if len(species) == 1: species = species[0]
                if atArrow:		products.append(species)
                else:			reactants.append(species)
                species = []
                if item == '->':  atArrow = True
            else:
                # Check that all reactant structures are in dictionary
                # The product structures are generated automatically and need
                # not be included
                if item not in self.dictionary and not atArrow:
                    raise InvalidDatabaseError('Reaction family template contains an unknown structure.')
                species.append(item)

        self.forwardTemplate = Reaction(reactants=reactants, products=products)

        # Second line is either "forward" (if reverse template different from
        # forward template) or "thermo_consistence" (if reverse template is the
        # same as forward template
        type = lines[1].strip().lower()
        index = 3
        if type == 'forward':
            index = 4

        # Remaining lines are reaction recipe for forward reaction
        self.forwardRecipe = ReactionRecipe()
        for line in lines[index:]:
            line = line.strip()

            # First item is the name of the action
            items = line.split()
            action = [ items[1].upper() ]
            assert action[0] in ['CHANGE_BOND','FORM_BOND','BREAK_BOND','GAIN_RADICAL','LOSE_RADICAL']

            # Remaining items are comma-delimited list of parameters enclosed by
            # {}, which we will split into individual parameters
            action.extend(items[2].strip()[1:-1].split(','))

            self.forwardRecipe.addAction(action)

        # Generate the reverse template
        if type == 'forward':
            self.generateProductTemplate()
            self.reverseTemplate = Reaction(reactants=self.forwardTemplate.products, products=self.forwardTemplate.reactants)
            self.reverseRecipe = self.forwardRecipe.getReverse()
        else:
            self.reverseTemplate = None
            self.reverseRecipe = None

    def generateProductTemplate(self):
        """
        Generate the product structures by applying the reaction template to
        the top-level nodes. For reactants defined by multiple structures, only
        the first is used here; it is assumed to be the most generic.
        """

        # First, generate a list of reactant structures that are actual
        # structures, rather than unions
        reactantStructures = []

        logging.log(0, "Generating template for products.")
        for reactant in self.forwardTemplate.reactants:
            if isinstance(reactant, list):	reactants = [reactant[0]]
            else:							reactants = [reactant]

            logging.log(0, "Reactants:%s"%reactants)
            for s in reactants: #
                struct = self.dictionary[s]
                #logging.debug("Reactant %s is class %s"%(str(s),struct.__class__))
                if isinstance(struct, LogicNode):
                    all_structures = struct.getPossibleStructures(self.dictionary)
                    logging.log(0, 'Expanding node %s to %s'%(s, all_structures))
                    reactantStructures.append(all_structures)
                else:
                    reactantStructures.append([struct])

        # Second, get all possible combinations of reactant structures
        reactantStructures = getAllCombinations(reactantStructures)

        # Third, generate all possible product structures by applying the
        # recipe to each combination of reactant structures
        # Note that bimolecular products are split by labeled atoms
        productStructures = []
        for reactantStructure in reactantStructures:
            productStructure = self.applyRecipe(reactantStructure, unique=False)
            productStructures.append(productStructure)

        # Fourth, remove duplicates from the lists
        productStructureList = [[] for i in range(len(productStructures[0]))]
        for productStructure in productStructures:
            for i, struct in enumerate(productStructure):
                for s in productStructureList[i]:
                    try:
                        if s.isIsomorphic(struct): break
                    except KeyError:
                        print struct.toAdjacencyList()
                        print s.toAdjacencyList()
                        raise
                else:
                    productStructureList[i].append(struct)
        
        # Fifth, associate structures with product template
        for i in range(len(self.forwardTemplate.products)):
            if len(productStructureList[i]) == 1:
                self.dictionary[self.forwardTemplate.products[i]] = productStructureList[i][0]
                self.tree.parent[self.forwardTemplate.products[i]] = None
                self.tree.children[self.forwardTemplate.products[i]] = []
            else:
                children = []
                for j in range(len(productStructureList[i])):
                    label = '%s_%i' % (self.forwardTemplate.products[i], j+1)
                    self.dictionary[label] = productStructureList[i][j]
                    children.append(label)
                    self.tree.parent[label] = self.forwardTemplate.products[i]
                    self.tree.children[label] = []

                self.tree.parent[self.forwardTemplate.products[i]] = None
                self.tree.children[self.forwardTemplate.products[i]] = children

                self.dictionary[self.forwardTemplate.products[i]] = LogicOr(children,invert=False)

    def reactantMatch(self, reactant, templateReactant):
        """
        Return ``True`` if the provided reactant matches the provided
        template reactant and ``False`` if not, along with a complete list of
        the identified mappings.
        """
        mapsList = []
        if templateReactant.__class__ == list: templateReactant = templateReactant[0]
        struct = self.dictionary[templateReactant]

        if isinstance(struct, LogicNode):
            for child_structure in struct.getPossibleStructures(self.dictionary):
                ismatch, mappings = reactant.findSubgraphIsomorphisms(child_structure)
                if ismatch:
                    mapsList.extend(mappings)
            return len(mapsList) > 0, mapsList
        elif isinstance(struct, Molecule):
            return reactant.findSubgraphIsomorphisms(struct)

    def applyRecipe(self, reactantStructures, forward=True, unique=True):
        """
        Apply the recipe for this reaction family to the list of
        :class:`Molecule` objects `reactantStructures`. The atoms
        of the reactant structures must already be tagged with the appropriate
        labels. Returns a list of structures corresponding to the products
        after checking that the correct number of products was produced.
        """

        # There is some hardcoding of reaction families in this function, so
        # we need the label of the reaction family for this
        label = self.label.lower()

        # Merge reactant structures into single structure
        # Also copy structures so we don't modify the originals
        # Since the tagging has already occurred, both the reactants and the
        # products will have tags
        if isinstance(reactantStructures[0], MoleculePattern):
            reactantStructure = MoleculePattern()
        else:
            reactantStructure = Molecule()
        for s in reactantStructures:
            reactantStructure = reactantStructure.merge(s.copy(deep=True))
        
        # Hardcoding of reaction family for radical recombination (colligation)
        # because the two reactants are identical, they have the same tags
        # In this case, we must change the labels from '*' and '*' to '*1' and
        # '*2'
        if label == 'r_recombination' and forward:
            identicalCenterCounter = 0
            for atom in reactantStructure.atoms:
                if atom.label == '*':
                    identicalCenterCounter += 1
                    atom.label = '*' + str(identicalCenterCounter)
            if identicalCenterCounter != 2:
                raise Exception('Unable to change labels from "*" to "*1" and "*2" for reaction family %s.' % (label))

        # Generate the product structure by applying the recipe
        if forward:
            self.forwardRecipe.applyForward(reactantStructure, unique)
        else:
            self.reverseRecipe.applyForward(reactantStructure, unique)
        productStructure = reactantStructure

        # Hardcoding of reaction family for reverse of radical recombination
        # (Unimolecular homolysis)
        # Because the two products are identical, they should the same tags
        # In this case, we must change the labels from '*1' and '*2' to '*' and
        # '*'
        if label == 'r_recombination' and not forward:
            for atom in productStructure.atoms:
                if atom.label == '*1' or atom.label == '*2': atom.label = '*'

        # If reaction family is its own reverse, relabel atoms
        if not self.reverseTemplate:
            # Get atom labels for products
            atomLabels = {}
            for atom in productStructure.atoms:
                if atom.label != '':
                    atomLabels[atom.label] = atom

            # This is hardcoding of reaction families (bad!)
            label = self.label.lower()
            if label == 'h_abstraction':
                # '*2' is the H that migrates
                # it moves from '*1' to '*3'
                atomLabels['*1'].label = '*3'
                atomLabels['*3'].label = '*1'

            elif label == 'intra_h_migration':
                # '*3' is the H that migrates
                # swap the two ends between which the H moves
                atomLabels['*1'].label = '*2'
                atomLabels['*2'].label = '*1'
                # reverse all the atoms in the chain between *1 and *2
                # i.e. swap *4 with the highest, *5 with the second-highest
                highest = len(atomLabels)
                if highest>4:
                    for i in range(4,highest+1):
                        atomLabels['*%d'%i].label = '*%d'%(4+highest-i)

        if not forward: template = self.reverseTemplate
        else:           template = self.forwardTemplate

        # Split product structure into multiple species if necessary
        if len(template.products) > 1:
            productStructures = productStructure.split()
        else:
            productStructures = [productStructure]

        for product in productStructures:
            product.updateConnectivityValues()

        # Make sure we've made the expected number of products
        if len(template.products) != len(productStructures):
            # We have a different number of products than expected by the template.
            # It might be because we found a ring-opening using a homolysis template
            if (label=='r_recombination' and not forward
             and len(productStructures) == 1
             and len(reactantStructures) == 1):
                # just be absolutely sure (maybe slow, but safe)
                rs = reactantStructures[0]
                if ( rs.isVertexInCycle(rs.getLabeledAtom('*1'))
                 and rs.isVertexInCycle(rs.getLabeledAtom('*2'))):
                    # both *1 and *2 are in cycles (probably the same one)
                    # so it's pretty safe to just fail quietly,
                    # and try the next reaction
                    return None

            # no other excuses, raise an exception
            message = 'Application of reaction recipe failed; expected %s product(s), but %s found.\n' % (len(self.template.products), len(productStructures))
            message += "Reaction family: %s \n"%str(self)
            message += "Reactant structures: %s \n"%reactantStructures
            message += "Product structures: %s \n"%productStructures
            message += "Template: %s"%self.template
            logging.error(message)
            #return None # don't fail!!! muhahaha
            raise Exception(message)

        # If there are two product structures, place the one containing '*1' first
        if len(productStructures) == 2:
            if not productStructures[0].containsLabeledAtom('*1') and \
                productStructures[1].containsLabeledAtom('*1'):
                productStructures.reverse()

        # If product structures are Molecule objects, update their atom types
        for struct in productStructures:
            if isinstance(struct, Molecule):
                struct.updateAtomTypes()

        # Return the product structures
        return productStructures

    def generateProductStructures(self, reactantStructures, maps):
        """
        For a given set of `reactantStructures` and a given set of `maps`,
        generate and return the corresponding product structures. The
        `reactantStructures` parameter should be given in the order the
        reactants are stored in the reaction family template. The `maps`
        parameter is a list of mappings of the top-level tree node of each
        *template* reactant to the corresponding *structure*. This function
        returns the product structures, species, and a boolean that tells
        whether any species are new.
        """

        # Clear any previous atom labeling from all reactant structures
        for struct in reactantStructures: struct.clearLabeledAtoms()

        # If there are two structures and they are the same, then make a copy
        # of the second one and adjust the second map to point to its atoms
        # This is for the case where A + A --> products
        if len(reactantStructures) == 2 and reactantStructures[0] == reactantStructures[1]:
            reactantStructures[1], newMap = reactantStructures[1].copy(returnMap=True)
            maps[1] = dict([(templateAtom,newMap[reactantAtom]) for templateAtom, reactantAtom in maps[1].iteritems()])

        # Tag atoms with labels
        for map in maps:
            for templateAtom, reactantAtom in map.iteritems():
                reactantAtom.label = templateAtom.label

        # Generate the product structures by applying the forward reaction recipe
        try:
            productStructures = self.applyRecipe(reactantStructures)
            if not productStructures: return None
        except chem.InvalidChemicalActionException, e:
            print 'Unable to apply reaction recipe!'
            print 'Reaction family is %s' % self
            print 'Reactant structures are:'
            for struct in reactantStructures:
                print struct.toAdjacencyList()
            raise

        # Check that reactant and product structures are allowed in this family
        # If not, then stop
        if self.forbidden is not None:
            for label, struct2 in self.forbidden.iteritems():
                for struct in reactantStructures:
                    if struct.isSubgraphIsomorphic(struct2): return None
                for struct in productStructures:
                    if struct.isSubgraphIsomorphic(struct2): return None

        return productStructures

    def createReaction(self, reactants, reactantStructures, productStructures, reactantAtomLabels):
        """
        Create a :class:`Reaction` object representing the reaction that
        converts the structures in `reactantStructures` corresponding to the
        species in `reactants` to the structures in `productStructures`. The
        atom labels for the reactants should already be known, and they are
        passed in the `reactantAtomLabels` parameter.
        """

        productAtomLabels = []
        for struct in productStructures:
            productAtomLabels.append(struct.getLabeledAtoms())

        # Convert product structures to product species
        products = []
        for i, struct0 in enumerate(productStructures):
            found, spec, struct, map = species.checkForExistingSpecies(struct0)
            if found:
                # Adjust atom labels mapping accordingly
                for label, atom in productAtomLabels[i].iteritems():
                    productAtomLabels[i][label] = map[atom]
                # Save struct rather than struct0
                productStructures[i] = struct
                # Append product species to list of products
                products.append(spec)
            else:
                product, isNew = species.makeNewSpecies(struct0, checkExisting=False)
                # Don't make a new reaction if no species was returned from
                # makeNewSpecies() (e.g. due to forbidden structure)
                if product is None: return None
                products.append(product)

        # Sort reactants and products (to make comparisons easier/faster)
        reactants.sort()
        products.sort()

        # Check that the reaction actually results in a different set of species
        if reactants == products:
            return None

        # Create reaction object
        forward = reaction.Reaction(reactants=reactants, products=products, family=self)
        reverse = reaction.Reaction(reactants=products, products=reactants, family=self.reverse or self)
        forward.reverse = reverse
        reverse.reverse = forward

        # Dictionaries containing the labeled atoms for the reactants and products
        forward.atomLabels = reactantAtomLabels
        reverse.atomLabels = productAtomLabels

        # Return the created reaction (forward direction only)
        return forward

    def getReactionList(self, reactants):
        """
        Generate a list of all of the possible reactions of this family between
        the list of `reactants`. The number of reactants provided must match
        the number of reactants expected by the template, or this function
        will return an empty list.
        """

        rxnList = []

        # Unimolecular reactants: A --> products
        if len(reactants) == 1 and self.template.isUnimolecular():

            # Iterate over all resonance isomers of the reactant
            for structure in reactants[0].structure:

                ismatch, map21, map12 = self.reactantMatch(structure, self.template.reactants[0])
                if ismatch:
                    for map in map12:

                        reactantAtomLabels = [{}]
                        for atom1, atom2 in map.iteritems():
                            reactantAtomLabels[0][atom1.label] = atom2

                        reactantStructures = [structure]
                        productStructures = self.generateProductStructures(reactantStructures, [map])
                        if productStructures:
                            rxn = self.createReaction(reactants, reactantStructures, productStructures, reactantAtomLabels)
                            if rxn: rxnList.append(rxn)

        # Bimolecular reactants: A + B --> products
        elif len(reactants) == 2 and self.template.isBimolecular():

            structuresA = reactants[0].structure
            structuresB = reactants[1].structure

            # Iterate over all resonance isomers of the reactant
            for structureA in structuresA:
                for structureB in structuresB:

                    # Reactants stored as A + B
                    ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[0])
                    ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[1])

                    # Iterate over each pair of matches (A, B)
                    if ismatch_A and ismatch_B:
                        for mapA in map12_A:
                            for mapB in map12_B:

                                reactantAtomLabels = [{},{}]
                                for atom1, atom2 in mapA.iteritems():
                                    reactantAtomLabels[0][atom1.label] = atom2
                                for atom1, atom2 in mapB.iteritems():
                                    reactantAtomLabels[1][atom1.label] = atom2

                                reactantStructures = [structureA, structureB]
                                productStructures = self.generateProductStructures(reactantStructures, [mapA, mapB])
                                if productStructures:
                                    rxn = self.createReaction(reactants, reactantStructures, productStructures, reactantAtomLabels)
                                    if rxn: rxnList.append(rxn)

                    # Only check for swapped reactants if they are different
                    if reactants[0].id != reactants[1].id:

                        # Reactants stored as B + A
                        ismatch_A, map21_A, map12_A = self.reactantMatch(structureA, self.template.reactants[1])
                        ismatch_B, map21_B, map12_B = self.reactantMatch(structureB, self.template.reactants[0])

                        # Iterate over each pair of matches (A, B)
                        if ismatch_A and ismatch_B:
                            for mapA in map12_A:
                                for mapB in map12_B:

                                    reactantAtomLabels = [{},{}]
                                    for atom1, atom2 in mapA.iteritems():
                                        reactantAtomLabels[0][atom1.label] = atom2
                                    for atom1, atom2 in mapB.iteritems():
                                        reactantAtomLabels[1][atom1.label] = atom2

                                    reactantStructures = [structureA, structureB]
                                    productStructures = self.generateProductStructures(reactantStructures, [mapA, mapB])
                                    if productStructures:
                                        rxn = self.createReaction(reactants, reactantStructures, productStructures, reactantAtomLabels)
                                        if rxn: rxnList.append(rxn)

        # Merge duplicate reactions and increment multiplier
        # In this context we already know that the family and the reactants
        # match, so we only need to check the products
        reactionsToRemove = []
        for i, rxn1 in enumerate(rxnList):
            for j, rxn2 in enumerate(rxnList[i+1:]):
                if rxn2 not in reactionsToRemove:
                    if rxn1.products == rxn2.products:
                        reactionsToRemove.append(rxn2)
                        rxn1.multiplier += 1.0
        for rxn in reactionsToRemove:
            rxnList.remove(rxn)

        # For R_Recombination reactions, the multiplier is twice what it should
        # be, so divide those by two
        # This is hardcoding of reaction families!
        if self.label.lower() == 'unimolecular homolysis':
            for rxn in rxnList:
                rxn.multiplier /= 2

        # Formally make the new reactions
        reactionsToRemove = []
        for i in range(len(rxnList)):
            rxn, isNew = reaction.makeNewReaction(rxnList[i])
            if isNew:
                rxnList[i] = rxn
            else:
                reactionsToRemove.append(rxnList[i])
        for rxn in reactionsToRemove:
            rxnList.remove(rxn)

        return rxnList

    def getKinetics(self, reaction, structures):
        """
        Determine the appropriate kinetics for `reaction` which involves the
        labeled atoms in `atoms`.
        """

        # Get forward reaction template and remove any duplicates
        forwardTemplate = self.tree.top[:]

        temporary=[]
        symmetric_tree=False
        for node in forwardTemplate:
            if node not in temporary:
                temporary.append(node)
            else:
                # duplicate node found at top of tree
                # eg. R_recombination: ['Y_rad', 'Y_rad']
                assert len(forwardTemplate)==2 , 'Can currently only do symmetric trees with nothing else in them'
                symmetric_tree = True
        forwardTemplate = temporary

        # Descend reactant trees as far as possible
        template = []
        for forward in forwardTemplate:
            # 'forward' is a head node that should be matched.
            # Get labeled atoms of forward
            node = forward
            group = self.dictionary[node]
            # to sort out "union" groups:
            # descends to the first child that's not a logical node
            if isinstance(group,LogicNode):
                all_structures = group.getPossibleStructures(self.dictionary)
                group = all_structures[0]
                node = 'First sub-structure of '+node
            # ...but this child may not match the structure.
            # eg. an R3 ring node will not match an R4 ring structure.
            # (but at least the first such child will contain fewest labels - we hope)

            atomList = group.getLabeledAtoms() # list of atom labels in highest non-union node

            for struct in structures:
                # Match labeled atoms
                # Check this structure has each of the atom labels in this group
                if not all([struct.containsLabeledAtom(label) for label in atomList]):
                    continue # don't try to match this structure - the atoms aren't there!
                # Match structures
                atoms = struct.getLabeledAtoms()
                matched_node = self.descendTree(struct, atoms, root=forward)
                if matched_node is not None:
                    template.append(matched_node)
                else:
                    logging.warning("Couldn't find match for %s in %s"%(forward,atomList))
                    logging.warning( struct.toAdjacencyList() )

        # Get fresh templates (with duplicate nodes back in)
        forwardTemplate = self.tree.top[:]
        if self.label.lower() == 'r_recombination':
            forwardTemplate.append(forwardTemplate[0])

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forwardTemplate is a list of the top level nodes that should be matched
        if len(template) != len(forwardTemplate):
            logging.warning('Unable to find matching template for reaction %s in reaction family %s' % (str(reaction), str(self)) )
            logging.warning(" Trying to match " + str(forwardTemplate))
            logging.warning(" Matched "+str(template))
            raise UndeterminableKineticsError(reaction)
            print str(self), template, forwardTemplate, reverseTemplate
            for reactant in reaction.reactants:
                print reactant.toAdjacencyList() + '\n'
            for product in reaction.products:
                print product.toAdjacencyList() + '\n'

            ## If unable to match template, use the most general template
            #template = forwardTemplate

        # climb the tree finding ancestors
        nodeLists = []
        for temp in template:
            nodeList = []
            while temp is not None:
                nodeList.append(temp)
                temp = self.tree.parent[temp]
            nodeLists.append(nodeList)

        # Generate all possible combinations of nodes
        items = getAllCombinations(nodeLists)

        # Generate list of kinetics at every node
        #logging.debug("   Template contains %s"%forwardTemplate)
        kinetics = []
        for item in items:
            itemData = self.library.getData(item)
            #logging.debug("   Looking for %s found %r"%(item, itemData))
            if itemData is not None:
                kinetics.append(itemData)

            if symmetric_tree: # we might only store kinetics the other way around
                item.reverse()
                itemData = self.library.getData(item)
                #logging.debug("   Also looking for %s found %r"%(item, itemData))
                if itemData is not None:
                    kinetics.append(itemData)

        # Make sure we've found at least one set of valid kinetics
        if len(kinetics) == 0:
            for reactant in structures:
                print reactant.toAdjacencyList() + '\n'
            raise UndeterminableKineticsError(reaction)

        # Choose the best kinetics
        # For now just return the kinetics with the highest index
        maxIndex = max([k.index for k in kinetics])
        kinetics = [k for k in kinetics if k.index == maxIndex][0]

        return kinetics.model

################################################################################

class KineticsGroupDatabase:
    """
    Represent a kinetics database, divided into reaction families. The
    `families` attribute stores a dictionary of :class:`ReactionFamily` objects
    representing the families in the set.
    """

    def __init__(self, path='', only_families=False):
        self.families = {}
        if path != '': self.load(path, only_families)

    def load(self, path, only_families=False):
        """
        Load a set of reaction families from the general database
        specified at `path`. If only_families is present, families not in
        this list will not be loaded (e.g. only_families=['H_Abstraction'] )
        """

        path = os.path.abspath(path)

        logging.info('Loading kinetics databases from %s...' % path)

        # Load the families from kinetics/families.txt
        familyList = []
        ffam = None
        try:
            ffam = open(os.path.join(path,'families.txt'), 'r')
            for line in ffam:
                line = removeCommentFromLine(line).strip()
                if len(line) > 0:
                    items = line.split()
                    items[0] = int(items[0])
                    familyList.append(items)
        except InvalidDatabaseError, e:
            logging.exception(str(e))
            return
        except IOError, e:
            logging.exception('Database file "' + e.filename + '" not found.')
            raise
        finally:
            if ffam: ffam.close()

        # Load the reaction families (if they exist and status is 'on')
        self.families = {}
        for index, status, label in familyList:
            fpath = os.path.join(path, label)
            if os.path.isdir(fpath) and status.lower() == 'on':
                # skip families not in only_families, if it's set
                if only_families and label not in only_families: continue

                logging.info('Loading reaction family %s...' % (label))
                family = ReactionFamily(label)
                family.load(fpath)
                self.families[family.label] = family
            else: logging.info("NOT loading reaction family %s." % label)

        logging.info('')

    def generateKineticsData(self, rxn, family, structures):
        return self.families[family].getKinetics(rxn, structures)

################################################################################

class KineticsPrimaryDatabase(Database):
    """
    A primary kinetics database, consisting of a dictionary of species
    and a library of reactions with corresponding kinetic data. (No tree is
    utilized in this database.) The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `label`             ``str``             A unique string identifier for the kinetics database
    `database`          :class:`Database`   Used to store a dictionary of species
    `reactions`         ``list``            A list of the reactions in the database
    `seedMechanism`     ``bool``            ``True`` to use as seed mechanism, ``False`` if not
    `reactionLibrary`   ``bool``            ``True`` to use as reaction library, ``False`` if not
    =================== =================== ====================================

    """

    def __init__(self, path=''):
        self.label = ''
        self.reactions = []
        self.seedMechanism = False
        self.reactionLibrary = False
        if path != '':
            self.load(path)
        else:
            self.database = None

    def __str__(self):
        return self.label

    def load(self, path, seedMechanism=False, reactionLibrary=False, old=False):
        """
        Load a primary kinetics database from the location `path`. Also set
        whether the kinetics database represents a `seedMechanism` or
        `reactionLibrary` (in addition to the default kinetics library).
        """
        self.label = os.path.split(path)[1]
        self.seedMechanism = seedMechanism
        self.reactionLibrary = reactionLibrary
        path = os.path.abspath(path)
        if seedMechanism:
            logging.info('Loading seed mechanism from %s...' % path)
        elif reactionLibrary:
            logging.info('Loading primary reaction database from %s...' % path)
        else:
            logging.info('Loading primary kinetics database from %s...' % path)
        if old:
            self.database = Database()
            self.__loadOldSpecies(os.path.join(path,'species.txt'))
            self.__loadOldReactions(os.path.join(path,'reactions.txt'))
        else:
            self.database = self.loadDatabase(os.path.join(path, 'database.py'))
        logging.info('')

    def __loadOldSpecies(self, path):
        """
        Load an old-style reaction library species file from `path`.
        """
        self.database.dictionary.load(path, pattern=False)
    
    def __loadOldReactions(self, path):
        """
        Load an old-style high-pressure limit reaction library from `path`.
        """
        self.reactions = []

        # Process the reactions file
        try:
            inUnitSection = False; inReactionSection = False
            Aunits = ''; Eunits = ''

            fdict = open(path, 'r')
            for line in fdict:
                line = removeCommentFromLine(line).strip()
                if len(line) > 0:
                    if inUnitSection:
                        if 'A:' in line or 'E:' in line:
                            units = line.split()[1]
                            if 'A:' in line:
                                Aunits = units.split('/') # Assume this is a 3-tuple: moles or molecules, volume, time
                                Aunits[1] = Aunits[1][0:-1] # Remove '3' from e.g. 'm3' or 'cm3'; this is assumed
                            elif 'E:' in line:
                                Eunits = units
                    elif inReactionSection:
                        reactants = []; products = []
                        items = line.split()
                        arrow = '=>' if '=>' in items else '<=>'
                        for item in items[0:items.index(arrow)]:
                            if item != '+':
                                if item not in self.database.dictionary:
                                    raise InvalidDatabaseError('Reactant %s not found in dictionary.' % item)
                                reactants.append(self.database.dictionary[item])
                        for item in items[items.index(arrow)+1:-6]:
                            if item != '+':
                                if item not in self.database.dictionary:
                                    raise InvalidDatabaseError('Product %s not found in dictionary.' % item)
                                products.append(self.database.dictionary[item])
                        A = pq.Quantity(float(items[-6]), '%s^%g/%s^%g/%s' %
                            (Aunits[1], 3*(len(reactants)-1), Aunits[0], len(reactants)-1, Aunits[2]))
                        n = pq.Quantity(float(items[-5]), '')
                        Ea = pq.Quantity(float(items[-4]), Eunits)
                        kin = ArrheniusModel(
                            A=float(A.simplified),
                            n=float(n.simplified),
                            Ea=float(Ea.simplified),
                            T0=1.0
                        )
                        rxn = Reaction(
                            reactants=reactants,
                            products=products,
                            kinetics=kin,
                            reversible= arrow =='<=>',
                        )
                        self.reactions.append(rxn)
                    if 'Unit:' in line:
                        inUnitSection = True; inReactionSection = False
                    elif 'Reactions:' in line:
                        inUnitSection = False; inReactionSection = True

        except (InvalidDatabaseError, InvalidAdjacencyListError), e:
            logging.exception(str(e))
            raise
        except IOError, e:
            logging.exception('Database dictionary file "' + e.filename + '" not found.')
            raise
        finally:
            if fdict: fdict.close()

    def generateKineticsData(self, reaction):
        """
        Determine the kinetics data for the given `reaction`, if it exists in
        the library. The reaction must be in the same direction as it is in
        the library in order to have the kinetics returned.
        """
        for libraryReaction in self.reactions:
            # To match, the reactions must have exactly the same numbers of
            # reactants and products
            if (len(libraryReaction.reactants) == len(reaction.reactants) and
              len(libraryReaction.products) == len(reaction.products)):
                # Try to pair up reactants
                reactants = reaction.reactants[:]
                for libraryReactant in libraryReaction.reactants:
                    for reactant in reactants:
                        if any([molecule.isIsomorphic(libraryReactant) for molecule in reactant.molecule]):
                            reactants.remove(reactant) # So we don't match the same reactant multiple times
                            break
                if len(reactants) != 0: continue
                # Try to pair up products
                products = reaction.products[:]
                for libraryProduct in libraryReaction.products:
                    for product in products:
                        if any([molecule.isIsomorphic(libraryProduct) for molecule in product.molecule]):
                            products.remove(product) # So we don't match the same product multiple times
                            break
                if len(products) != 0: continue
                # If we're here, then we were able to match all reactants and
                # products, so we assume that the reactions match
                return libraryReaction.kinetics

        # If we're here, then there were no matches found, so return no match
        return None

    def getReactionList(self, reactants):
        """
        Generate a list of all of the possible reactions of this family between
        the list of `reactants`.
        """
        rxnList = []

        for libraryReaction in self.reactions:
            # To match, the reactions must have exactly the same numbers of
            # reactants
            if len(libraryReaction.reactants) == len(reactants):
                # Try to pair up reactants
                reactants0 = reactants[:]
                for libraryReactant in libraryReaction.reactants:
                    for reactant in reactants0:
                        if any([molecule.isIsomorphic(libraryReactant) for molecule in reactant.molecule]):
                            reactants0.remove(reactant) # So we don't match the same reactant multiple times
                            break
                if len(reactants0) == 0:
                    rxnList.append(libraryReaction)

            # To match, the reactions must have exactly the same numbers of
            # reactants
            if len(libraryReaction.products) == len(reactants):
                # Try to pair up products
                reactants0 = reactants[:]
                for libraryProduct in libraryReaction.products:
                    for reactant in reactants0:
                        if any([molecule.isIsomorphic(libraryProduct) for molecule in reactant.molecule]):
                            reactants0.remove(reactant) # So we don't match the same product multiple times
                            break
                if len(reactants0) == 0:
                    rxnList.append(libraryReaction)

        return rxnList

################################################################################

kineticsDatabases = []

def loadKineticsDatabase(dstr, group=True, old=False, seedMechanism=False, reactionLibrary=False, only_families=False):
    """
    Load the RMG kinetics database located at `dstr` into the global variable
    `rmg.reaction.kineticsDatabase`.
    """
    global kineticsDatabases
    if group:
        kineticsDatabase = KineticsGroupDatabase()
        kineticsDatabase.load(dstr, only_families)
        kineticsDatabases.append(kineticsDatabase)
    else:
        kineticsDatabase = KineticsPrimaryDatabase()
        kineticsDatabase.load(dstr, old=old, seedMechanism=seedMechanism, reactionLibrary=reactionLibrary)
        kineticsDatabases.append(kineticsDatabase)
    return kineticsDatabase

def generateKineticsData(rxn, family, structures):
    """
    Get the kinetics data associated with `reaction` in `family` by looking in
    the loaded kinetics database. The reactants and products in the molecule
    should already be labeled appropriately. An
    :class:`UndeterminableKineticsError` is raised if no kinetics could be
    determined for the given reaction.
    """
    global kineticsDatabases
    kinetics = None
    for kineticsDatabase in kineticsDatabases:
        if isinstance(kineticsDatabase, KineticsGroupDatabase):
            kinetics = kineticsDatabase.generateKineticsData(rxn, family, structures)
        elif isinstance(kineticsDatabase, KineticsPrimaryDatabase) and not kineticsDatabase.seedMechanism and not kineticsDatabase.reactionLibrary:
            # Assume that, for seed mechanisms and reaction libraries, you've
            # already determined the kinetics (if possible) when you discovered
            # the reaction
            kinetics = kineticsDatabase.generateKineticsData(rxn)
        if kinetics is not None: return kinetics
    raise UndeterminableKineticsError(rxn)
    
################################################################################

