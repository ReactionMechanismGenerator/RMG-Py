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

"""
This module contains functionality for working with kinetics libraries.
"""

import os
import os.path
import re
import logging
import codecs
from copy import copy, deepcopy

from rmgpy.data.base import *

from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.kinetics import Arrhenius, ArrheniusEP, ThirdBody, Lindemann, Troe, \
                           PDepArrhenius, MultiArrhenius, MultiPDepArrhenius, \
                           Chebyshev, KineticsData, PDepKineticsModel
from rmgpy.molecule import Bond, GroupBond, Group
from rmgpy.species import Species

from .common import KineticsError, saveEntry
from .depository import KineticsDepository

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

class ReactionPairsError(Exception):
    """
    An exception to be raised when an error occurs while working with reaction
    pairs.
    """
    pass

################################################################################

class TemplateReaction(Reaction):
    """
    A Reaction object generated from a reaction family template. In addition to
    the usual attributes, this class includes a `family` attribute to store the
    family that it was created from, as well as a `estimator` attribute to indicate
    whether it came from a rate rules or a group additivity estimate.
    """

    def __init__(self,
                index=-1,
                reactants=None,
                products=None,
                kinetics=None,
                reversible=True,
                transitionState=None,
                duplicate=False,
                degeneracy=1,
                pairs=None,
                family=None,
                template=None,
                estimator=None,
                ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          kinetics=kinetics,
                          reversible=reversible,
                          transitionState=transitionState,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs
                          )
        self.family = family
        self.template = template
        self.estimator = estimator

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TemplateReaction, (self.index,
                                   self.reactants,
                                   self.products,
                                   self.kinetics,
                                   self.reversible,
                                   self.transitionState,
                                   self.duplicate,
                                   self.degeneracy,
                                   self.pairs,
                                   self.family,
                                   self.template,
                                   self.estimator
                                   ))

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        TemplateReaction this should be a KineticsGroups object.
        """
        return self.family

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

        pattern = isinstance(struct, Group)

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
                    bond = GroupBond(atom1, atom2, order=['S']) if pattern else Bond(atom1, atom2, order='S')
                    struct.addBond(bond)
                    atom1.applyAction(['FORM_BOND', label1, info, label2])
                    atom2.applyAction(['FORM_BOND', label1, info, label2])
                elif (action[0] == 'BREAK_BOND' and doForward) or (action[0] == 'FORM_BOND' and not doForward):
                    if not struct.hasBond(atom1, atom2):
                        raise InvalidActionError('Attempted to remove a nonexistent bond.')
                    bond = struct.getBond(atom1, atom2)
                    struct.removeBond(bond)
                    atom1.applyAction(['BREAK_BOND', label1, info, label2])
                    atom2.applyAction(['BREAK_BOND', label1, info, label2])

            elif action[0] in ['LOSE_RADICAL', 'GAIN_RADICAL']:

                label, change = action[1:]
                change = int(change)

                # Find associated atom
                atom = struct.getLabeledAtom(label)
                if atom is None:
                    raise InvalidActionError('Unable to find atom with label "{0}" while applying reaction recipe.'.format(label))

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

class KineticsGroups(Database):
    """
    A class for working with an RMG kinetics family group additivity values. 
    """

    def __init__(self,
                 entries=None,
                 top=None,
                 label='',
                 name='',
                 shortDesc='',
                 longDesc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None
                 ):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
        self.numReactants = 0
        
    def __repr__(self):
        return '<KineticsGroups "{0}">'.format(self.label)

    def loadEntry(self, index, label, group, kinetics, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        if group[0:3].upper() == 'OR{' or group[0:4].upper() == 'AND{' or group[0:7].upper() == 'NOT OR{' or group[0:8].upper() == 'NOT AND{':
            item = makeLogicNode(group)
        else:
            item = Group().fromAdjacencyList(group)
        self.entries[label] = Entry(
            index = index,
            label = label,
            item = item,
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )

    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """

        # Get forward reaction template and remove any duplicates
        forwardTemplate = self.top[:]

        temporary = []
        symmetricTree = False
        for entry in forwardTemplate:
            if entry not in temporary:
                temporary.append(entry)
            else:
                # duplicate node found at top of tree
                # eg. R_recombination: ['Y_rad', 'Y_rad']
                assert len(forwardTemplate)==2 , 'Can currently only do symmetric trees with nothing else in them'
                symmetricTree = True
        forwardTemplate = temporary

        # Descend reactant trees as far as possible
        template = []
        for entry in forwardTemplate:
            # entry is a top-level node that should be matched
            group = entry.item

            # To sort out "union" groups, descend to the first child that's not a logical node
            # ...but this child may not match the structure.
            # eg. an R3 ring node will not match an R4 ring structure.
            # (but at least the first such child will contain fewest labels - we hope)
            if isinstance(entry.item, LogicNode):
                group = entry.item.getPossibleStructures(self.entries)[0]

            atomList = group.getLabeledAtoms() # list of atom labels in highest non-union node

            for reactant in reaction.reactants:
                if isinstance(reactant, Species):
                    reactant = reactant.molecule[0]
                # Match labeled atoms
                # Check this reactant has each of the atom labels in this group
                if not all([reactant.containsLabeledAtom(label) for label in atomList]):
                    continue # don't try to match this structure - the atoms aren't there!
                # Match structures
                atoms = reactant.getLabeledAtoms()
                matched_node = self.descendTree(reactant, atoms, root=entry)
                if matched_node is not None:
                    template.append(matched_node)
                #else:
                #    logging.warning("Couldn't find match for {0} in {1}".format(entry,atomList))
                #    logging.warning(reactant.toAdjacencyList())

        # Get fresh templates (with duplicate nodes back in)
        forwardTemplate = self.top[:]
        if self.label.lower().startswith('r_recombination'):
            forwardTemplate.append(forwardTemplate[0])

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forwardTemplate is a list of the top level nodes that should be matched
        if len(template) != len(forwardTemplate):
            #logging.warning('Unable to find matching template for reaction {0} in reaction family {1}'.format(str(reaction), str(self)) )
            #logging.warning(" Trying to match " + str(forwardTemplate))
            #logging.warning(" Matched "+str(template))
            #print str(self), template, forwardTemplate
            #for reactant in reaction.reactants:
            #    print reactant.toAdjacencyList() + '\n'
            #for product in reaction.products:
            #    print product.toAdjacencyList() + '\n'
            raise UndeterminableKineticsError(reaction)

        return template

    def estimateKineticsUsingGroupAdditivity(self, template, referenceKinetics, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using group additivity.
        """

        # Start with the generic kinetics of the top-level nodes
        # Make a copy so we don't modify the original
        kinetics = deepcopy(referenceKinetics)
        
        # Now add in more specific corrections if possible
        for node in template:
            entry = node
            comment_line = "Matched node "
            while entry.data is None and entry not in self.top:
                # Keep climbing tree until you find a (non-top) node with data.
                comment_line += "{0} >> ".format(entry.label)
                entry = entry.parent
            if entry.data is not None and entry not in self.top:
                kinetics = self.__multiplyKineticsData(kinetics, entry.data)
                comment_line += "{0} ({1})".format(entry.label, entry.longDesc.split('\n')[0])
            elif entry in self.top:
                comment_line += "{0} (Top node)".format(entry.label)
            kinetics.comment += comment_line + '\n'

        # Also include reaction-path degeneracy
        if isinstance(kinetics, KineticsData):
            kinetics.kdata.value_si *= degeneracy
        elif isinstance(kinetics, Arrhenius):
            kinetics.A.value_si *= degeneracy
        elif kinetics is not None:
            raise KineticsError('Unexpected kinetics type "{0}" encountered while generating kinetics from group values.'.format(kinetics.__class__))
        kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)
        
        return kinetics

    def __multiplyKineticsData(self, kinetics1, kinetics2):
        """
        Multiply two kinetics objects `kinetics1` and `kinetics2` of the same
        class together, returning their product as a new kinetics object of 
        that class. Currently this only works for :class:`KineticsData` or
        :class:`Arrhenius` objects.
        """
        if isinstance(kinetics1, KineticsData) and isinstance(kinetics2, KineticsData):
            if len(kinetics1.Tdata.value_si) != len(kinetics2.Tdata.value_si) or any([T1 != T2 for T1, T2 in zip(kinetics1.Tdata.value_si, kinetics2.Tdata.value_si)]):
                raise KineticsError('Cannot add these KineticsData objects due to their having different temperature points.')
            kinetics = KineticsData(
                Tdata = (kinetics1.Tdata.value, kinetics2.Tdata.units),
                kdata = (kinetics1.kdata.value * kinetics2.kdata.value, kinetics1.kdata.units),
            )
        elif isinstance(kinetics1, Arrhenius) and isinstance(kinetics2, Arrhenius):
            assert kinetics1.A.units == kinetics2.A.units
            assert kinetics1.Ea.units == kinetics2.Ea.units
            assert kinetics1.T0.units == kinetics2.T0.units
            assert kinetics1.T0.value == kinetics2.T0.value
            kinetics = Arrhenius(
                A = (kinetics1.A.value * kinetics2.A.value, kinetics1.A.units),
                n = (kinetics1.n.value + kinetics2.n.value, kinetics1.n.units),
                Ea = (kinetics1.Ea.value + kinetics2.Ea.value, kinetics1.Ea.units),
                T0 = (kinetics1.T0.value, kinetics1.T0.units),
            )
        else:
            raise KineticsError('Unable to multiply kinetics types "{0}" and "{1}".'.format(kinetics1.__class__, kinetics2.__class__))
        
        if kinetics1.Tmin is not None and kinetics2.Tmin is not None:
            kinetics.Tmin = kinetics1.Tmin if kinetics1.Tmin.value_si > kinetics2.Tmin.value_si else kinetics2.Tmin
        elif kinetics1.Tmin is not None and kinetics2.Tmin is None:
            kinetics.Tmin = kinetics1.Tmin
        elif kinetics1.Tmin is None and kinetics2.Tmin is not None:
            kinetics.Tmin = kinetics2.Tmin
        
        if kinetics1.Tmax is not None and kinetics2.Tmax is not None:
            kinetics.Tmax = kinetics1.Tmax if kinetics1.Tmax.value_si < kinetics2.Tmax.value_si else kinetics2.Tmax
        elif kinetics1.Tmax is not None and kinetics2.Tmax is None:
            kinetics.Tmax = kinetics1.Tmax
        elif kinetics1.Tmax is None and kinetics2.Tmax is not None:
            kinetics.Tmax = kinetics2.Tmax
        
        if kinetics1.Pmin is not None and kinetics2.Pmin is not None:
            kinetics.Pmin = kinetics1.Pmin if kinetics1.Pmin.value_si > kinetics2.Pmin.value_si else kinetics2.Pmin
        elif kinetics1.Pmin is not None and kinetics2.Pmin is None:
            kinetics.Pmin = kinetics1.Pmin
        elif kinetics1.Pmin is None and kinetics2.Pmin is not None:
            kinetics.Pmin = kinetics2.Pmin
        
        if kinetics1.Pmax is not None and kinetics2.Pmax is not None:
            kinetics.Pmax = kinetics1.Pmax if kinetics1.Pmax.value_si < kinetics2.Pmax.value_si else kinetics2.Pmax
        elif kinetics1.Pmax is not None and kinetics2.Pmax is None:
            kinetics.Pmax = kinetics1.Pmax
        elif kinetics1.Pmax is None and kinetics2.Pmax is not None:
            kinetics.Pmax = kinetics2.Pmax
        
        if kinetics1.comment == '': kinetics.comment = kinetics2.comment
        elif kinetics2.comment == '': kinetics.comment = kinetics1.comment
        else: kinetics.comment = kinetics1.comment + ' + ' + kinetics2.comment
        return kinetics

################################################################################

class KineticsFamily(Database):
    """
    A class for working with an RMG kinetics family: a set of reactions with 
    similar chemistry, and therefore similar reaction rates. The attributes 
    are:

    =================== =============================== ========================
    Attribute           Type                            Description
    =================== =============================== ========================
    `reverse`           ``string``                      The name of the reverse reaction family
    `forwardTemplate`   :class:`Reaction`               The forward reaction template
    `forwardRecipe`     :class:`ReactionRecipe`         The steps to take when applying the forward reaction to a set of reactants
    `reverseTemplate`   :class:`Reaction`               The reverse reaction template
    `reverseRecipe`     :class:`ReactionRecipe`         The steps to take when applying the reverse reaction to a set of reactants
    `forbidden`         :class:`ForbiddenStructures`    (Optional) Forbidden product structures in either direction
    `ownReverse`        `Boolean`                       It's its own reverse?
    ------------------- ------------------------------- ------------------------
    `groups`            :class:`KineticsGroups`         The set of kinetics group additivity values
    `rules`             :class:`KineticsDepository`     The depository of kinetics rate rules from RMG-Java
    `depositories`      ``dict``                        A set of additional depositories used to store kinetics data from various sources
    =================== =============================== ========================

    There are a few reaction families that are their own reverse (hydrogen
    abstraction and intramolecular hydrogen migration); for these
    `reverseTemplate` and `reverseRecipe` will both be ``None``.
    """

    def __init__(self,
                 entries=None,
                 top=None,
                 label='',
                 name='',
                 reverse='',
                 shortDesc='',
                 longDesc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None
                 ):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
        self.reverse = reverse
        self.forwardTemplate = forwardTemplate
        self.forwardRecipe = forwardRecipe
        self.reverseTemplate = reverseTemplate
        self.reverseRecipe = reverseRecipe
        self.forbidden = forbidden
        self.ownReverse = forwardTemplate is not None and reverseTemplate is None
        # Kinetics depositories of training and test data
        self.groups = None
        self.rules = None
        self.depositories = []

    def __repr__(self):
        return '<ReactionFamily "{0}">'.format(self.label)

    def loadOld(self, path):
        """
        Load an old-style RMG kinetics group additivity database from the
        location `path`.
        """
        self.label = os.path.basename(path)
        self.name = self.label

        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        self.groups.name = self.groups.label
        try:
            self.groups.loadOldDictionary(os.path.join(path, 'dictionary.txt'), pattern=True)
        except Exception:
            logging.error('Error while reading old kinetics family dictionary from {0!r}.'.format(path))
            raise
        try:
            self.groups.loadOldTree(os.path.join(path, 'tree.txt'))
        except Exception:
            logging.error('Error while reading old kinetics family tree from {0!r}.'.format(path))
            raise

        # The old kinetics groups use rate rules (not group additivity values),
        # so we can't load the old rateLibrary.txt
        
        # Load the reaction recipe
        try:
            self.loadOldTemplate(os.path.join(path, 'reactionAdjList.txt'))
        except Exception:
            logging.error('Error while reading old kinetics family template/recipe from {0!r}.'.format(path))
            raise
        # Construct the forward and reverse templates
        reactants = [self.groups.entries[label] for label in self.forwardTemplate.reactants]
        if self.ownReverse:
            self.forwardTemplate = Reaction(reactants=reactants, products=reactants)
            self.reverseTemplate = None
        else:
            products = self.generateProductTemplate(reactants)
            self.forwardTemplate = Reaction(reactants=reactants, products=products)
            self.reverseTemplate = Reaction(reactants=reactants, products=products)

        self.groups.numReactants = len(self.forwardTemplate.reactants)

        # Load forbidden structures if present
        try:
            if os.path.exists(os.path.join(path, 'forbiddenGroups.txt')):
                self.forbidden = ForbiddenStructures().loadOld(os.path.join(path, 'forbiddenGroups.txt'))
        except Exception:
            logging.error('Error while reading old kinetics family forbidden groups from {0!r}.'.format(path))
            raise
            
        entries = self.groups.top[:]
        for entry in self.groups.top:
            entries.extend(self.groups.descendants(entry))
        for index, entry in enumerate(entries):
            entry.index = index + 1
            
        self.rules = KineticsDepository(label='{0}/rules'.format(self.label),
                                        recommended=True)
        self.rules.name = self.rules.label
        try:
            self.rules.loadOldRateRules(path, self.groups, numLabels=max(len(self.forwardTemplate.reactants), len(self.groups.top)))
        except Exception:
            logging.error('Error while reading old kinetics family rules from {0!r}.'.format(path))
            raise
        self.depositories = {}

        return self

    def loadOldTemplate(self, path):
        """
        Load an old-style RMG reaction family template from the location `path`.
        """

        self.forwardTemplate = Reaction(reactants=[], products=[])
        self.forwardRecipe = ReactionRecipe()
        self.ownReverse = False

        ftemp = None
        # Process the template file
        try:
            ftemp = open(path, 'r')
            for line in ftemp:
                line = line.strip()
                if len(line) > 0 and line[0] == '(':
                    # This is a recipe action line
                    tokens = line.split()
                    action = [tokens[1]]
                    action.extend(tokens[2][1:-1].split(','))
                    self.forwardRecipe.addAction(action)
                elif 'thermo_consistence' in line:
                    self.ownReverse = True
                elif 'reverse' in line:
                    self.reverse = line.split(':')[1].strip()
                elif '->' in line:
                    # This is the template line
                    tokens = line.split()
                    atArrow = False
                    for token in tokens:
                        if token == '->':
                            atArrow = True
                        elif token != '+' and not atArrow:
                            self.forwardTemplate.reactants.append(token)
                        elif token != '+' and atArrow:
                            self.forwardTemplate.products.append(token)
        except IOError, e:
            logging.exception('Database template file "' + e.filename + '" not found.')
            raise
        finally:
            if ftemp: ftemp.close()

    def saveOld(self, path):
        """
        Save the old RMG kinetics groups to the given `path` on disk.
        """
        if not os.path.exists(path): os.mkdir(path)
        
        self.groups.saveOldDictionary(os.path.join(path, 'dictionary.txt'))
        self.groups.saveOldTree(os.path.join(path, 'tree.txt'))
        # The old kinetics groups use rate rules (not group additivity values),
        # so we can't save the old rateLibrary.txt
        self.saveOldTemplate(os.path.join(path, 'reactionAdjList.txt'))
        # Save forbidden structures if present
        if self.forbidden is not None:
            self.forbidden.saveOld(os.path.join(path, 'forbiddenGroups.txt'))
            
        self.rules.saveOldRateRules(path, self)
            
    def saveOldTemplate(self, path):
        """
        Save an old-style RMG reaction family template from the location `path`.
        """
        ftemp = open(path, 'w')
        
        # Write the template
        ftemp.write('{0} -> {1}\n'.format(
            ' + '.join([entry.label for entry in self.forwardTemplate.reactants]),
            ' + '.join([entry.label for entry in self.forwardTemplate.products]),
        ))
        ftemp.write('\n')
        
        # Write the reaction type and reverse name
        if self.ownReverse:
            ftemp.write('thermo_consistence\n')
        else:
            ftemp.write('forward\n')
            ftemp.write('reverse: {0}\n'.format(self.reverse))
        ftemp.write('\n')
        
        # Write the reaction recipe
        ftemp.write('Actions 1\n')
        for index, action in enumerate(self.forwardRecipe.actions):
            ftemp.write('({0}) {1:<15} {{{2}}}\n'.format(index+1, action[0], ','.join(action[1:])))
        ftemp.write('\n')
        
        ftemp.close()
    
    def load(self, path, local_context=None, global_context=None, depositoryLabels=None):
        """
        Load a kinetics database from a file located at `path` on disk.
        
        If `depositoryLabels` is a list, eg. ['training','PrIMe'], then only those
        depositories are loaded, and they are searched in that order when
        generating kinetics.
        
        If depositoryLabels is None then load 'training' first then everything else.
        If depositoryLabels is not None then load in the order specified in depositoryLabels.
        """
        local_context['recipe'] = self.loadRecipe
        local_context['template'] = self.loadTemplate
        local_context['forbidden'] = self.loadForbidden
        local_context['True'] = True
        local_context['False'] = False
        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        logging.debug("Loading kinetics family groups from {0}".format(os.path.join(path, 'groups.py')))
        Database.load(self.groups, os.path.join(path, 'groups.py'), local_context, global_context)
        self.name = self.label
        
        # Generate the reverse template if necessary
        self.forwardTemplate.reactants = [self.groups.entries[label] for label in self.forwardTemplate.reactants]
        if self.ownReverse:
            self.forwardTemplate.products = self.forwardTemplate.reactants[:]
            self.reverseTemplate = None
            self.reverseRecipe = None
        else:
            try:
                self.reverse = local_context['reverse']
            except KeyError:
                self.reverse = '{0}_reverse'.format(self.label)
            self.forwardTemplate.products = self.generateProductTemplate(self.forwardTemplate.reactants)
            self.reverseTemplate = Reaction(reactants=self.forwardTemplate.products, products=self.forwardTemplate.reactants)
            self.reverseRecipe = self.forwardRecipe.getReverse()
        
        self.groups.numReactants = len(self.forwardTemplate.reactants)
            
        self.rules = KineticsDepository(label='{0}/rules'.format(self.label))
        logging.debug("Loading kinetics family rules from {0}".format(os.path.join(path, 'rules.py')))
        self.rules.load(os.path.join(path, 'rules.py'), local_context, global_context)
        
        self.depositories = []
        # If depositoryLabels is None then load 'training' first then everything else.
        # If depositoryLabels is not None then load in the order specified in depositoryLabels.
        for name in (['training'] if depositoryLabels is None else depositoryLabels) :
            label = '{0}/{1}'.format(self.label, name)
            f = name+'.py'
            fpath = os.path.join(path,f)
            if not os.path.exists(fpath):
                logging.warning("Requested depository {0} does not exist".format(fpath))
                continue
            depository = KineticsDepository(label=label)
            logging.debug("Loading kinetics family depository from {0}".format(fpath))
            depository.load(fpath, local_context, global_context)
            self.depositories.append(depository)
        
        if depositoryLabels is None:
            # load all the remaining depositories, in order returned by os.walk
            for root, dirs, files in os.walk(path):
                if 'training' in root: continue
                for f in files:
                    if not f.endswith('.py'): continue
                    name = f.split('.py')[0]
                    if name not in ['groups', 'rules'] and name not in (depositoryLabels or ['training']):
                        fpath = os.path.join(root, f)
                        label = '{0}/{1}'.format(self.label, name)
                        depository = KineticsDepository(label=label)
                        logging.debug("Loading kinetics family depository from {0}".format(fpath))
                        depository.load(fpath, local_context, global_context)
                        self.depositories.append(depository)
            
    def loadTemplate(self, reactants, products, ownReverse=False):
        """
        Load information about the reaction template.
        """
        self.forwardTemplate = Reaction(reactants=reactants, products=products)
        self.ownReverse = ownReverse

    def loadRecipe(self, actions):
        """
        Load information about the reaction recipe.
        """
        # Remaining lines are reaction recipe for forward reaction
        self.forwardRecipe = ReactionRecipe()
        for action in actions:
            action[0] = action[0].upper()
            assert action[0] in ['CHANGE_BOND','FORM_BOND','BREAK_BOND','GAIN_RADICAL','LOSE_RADICAL']
            self.forwardRecipe.addAction(action)

    def loadForbidden(self, label, group, shortDesc='', longDesc='', history=None):
        """
        Load information about a forbidden structure.
        """
        if not self.forbidden:
            self.forbidden = ForbiddenStructures()
        self.forbidden.loadEntry(label=label, group=group, shortDesc=shortDesc, longDesc=longDesc, history=history)

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def save(self, path, entryName='entry'):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """
        self.saveGroups(os.path.join(path, 'groups.py'), entryName=entryName)
        self.rules.save(os.path.join(path, 'rules.py'))
        for label, depository in self.depositories.iteritems():
            self.saveDepository(depository, os.path.join(path, '{0}.py'.format(label[len(self.label)+1:])))
    
    def saveDepository(self, depository, path):
        """
        Save the given kinetics family `depository` to the location `path` on
        disk.
        """
        depository.save(os.path.join(path))        
        
    def saveGroups(self, path, entryName='entry'):
        """
        Save the current database to the file at location `path` on disk. The
        optional `entryName` parameter specifies the identifier used for each
        data entry.
        """
        entries = self.groups.getEntriesToSave()
                
        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}/groups"\n'.format(self.name))
        f.write('shortDesc = u"{0}"\n'.format(self.shortDesc))
        f.write('longDesc = u"""\n')
        f.write(self.longDesc)
        f.write('\n"""\n\n')

        # Write the template
        f.write('template(reactants=[{0}], products=[{1}], ownReverse={2})\n\n'.format(
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.reactants]),
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.products]),
            self.ownReverse))

        # Write reverse name
        if not self.ownReverse:
            f.write('reverse = "{0}"\n\n'.format(self.reverse))

        # Write the recipe
        f.write('recipe(actions=[\n')
        for action in self.forwardRecipe.actions:
            f.write('    {0!r},\n'.format(action))
        f.write('])\n\n')

        # Save the entries
        for entry in entries:
            self.saveEntry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generateOldTree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        # Save forbidden structures, if present
        if self.forbidden is not None:
            entries = self.forbidden.entries.values()
            entries.sort(key=lambda x: x.label)
            for entry in entries:
                self.forbidden.saveEntry(f, entry, name='forbidden')
    
        f.close()

    def generateProductTemplate(self, reactants0):
        """
        Generate the product structures by applying the reaction template to
        the top-level nodes. For reactants defined by multiple structures, only
        the first is used here; it is assumed to be the most generic.
        """

        # First, generate a list of reactant structures that are actual
        # structures, rather than unions
        reactantStructures = []

        logging.log(0, "Generating template for products.")
        for reactant in reactants0:
            if isinstance(reactant, list):  reactants = [reactant[0]]
            else:                           reactants = [reactant]

            logging.log(0, "Reactants: {0}".format(reactants))
            for s in reactants: #
                struct = s.item
                if isinstance(struct, LogicNode):
                    all_structures = struct.getPossibleStructures(self.groups.entries)
                    logging.log(0, 'Expanding node {0} to {1}'.format(s, all_structures))
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
            productStructure = self.applyRecipe(reactantStructure, forward=True, unique=False)
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
        productSet = []
        for index, products in enumerate(productStructureList):
            label = self.forwardTemplate.products[index]
            if len(products) == 1:
                entry = Entry(
                    label = label,
                    item = products[0],
                )
                self.groups.entries[entry.label] = entry
                productSet.append(entry)
            else:
                item = []
                counter = 0
                for product in products:
                    entry = Entry(
                        label = '{0}{1:d}'.format(label,counter+1),
                        item = product,
                    )
                    item.append(entry.label)
                    self.groups.entries[entry.label] = entry
                    counter += 1

                item = LogicOr(item,invert=False)
                entry = Entry(
                    label = label,
                    item = item,
                )
                self.groups.entries[entry.label] = entry
                counter += 1
                productSet.append(entry)

        return productSet

    def hasRateRule(self, template):
        """
        Return ``True`` if a rate rule with the given `template` currently 
        exists, or ``False`` otherwise.
        """
        try:
            return self.getRateRule(template) is not None
        except ValueError:
            return False

    def getRateRule(self, template):
        """
        Return the rate rule with the given `template`. Raises a 
        :class:`ValueError` if no corresponding entry exists.
        """
        entries = []
        templateLabels = ';'.join([group.label for group in template])
        for entry in self.rules.entries.values():
            if templateLabels == entry.label:
                entries.append(entry)
        
        if self.label.lower() == 'r_recombination' and template[0] != template[1]:
            template.reverse()
            templateLabels = ';'.join([group.label for group in template])
            for entry in self.rules.entries.values():
                if templateLabels == entry.label:
                    entries.append(entry)
            template.reverse()
            
        if len(entries) == 1:
            return entries[0]
        elif len(entries) > 1:
            if any([entry.rank > 0 for entry in entries]):
                entries = [entry for entry in entries if entry.rank > 0]
                entries.sort(key=lambda x: (x.rank, x.index))
                return entries[0]
            else:
                entries.sort(key=lambda x: x.index)
                return entries[0]
        else:
            raise ValueError('No entry for template {0}.'.format(template))

    def addKineticsRulesFromTrainingSet(self, thermoDatabase=None):
        """
        For each reaction involving real reactants and products in the training
        set, add a rate rule for that reaction.
        """
        for depository in self.depositories:
            if depository.label.endswith('training'):
                break
        else:
            raise Exception('Could not find training depository in family {0}.'.format(self.label))
        
        index = max([e.index for e in self.rules.entries.values()] or [0]) + 1
        
        entries = depository.entries.values()
        entries.sort(key=lambda x: x.index)
        reverse_entries = []
        for entry in entries:
            try:
                template = self.getReactionTemplate(entry.item)
            except UndeterminableKineticsError:
                # Some entries might be stored in the reverse direction for
                # this family; save them so we can try this
                reverse_entries.append(entry)
                continue
            
            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            
            new_entry = Entry(
                index = index,
                label = ';'.join([g.label for g in template]),
                item = template,
                data = ArrheniusEP(
                    A = deepcopy(data.A),
                    n = deepcopy(data.n),
                    alpha = 0,
                    E0 = deepcopy(data.Ea),
                    Tmin = deepcopy(data.Tmin),
                    Tmax = deepcopy(data.Tmax),
                ),
                rank = 3,
            )
            new_entry.data.A.value_si /= entry.item.degeneracy
            self.rules.entries[index] = new_entry
            index += 1
        
        # Process the entries that are stored in the reverse direction of the
        # family definition
        for entry in reverse_entries:
            
            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            
            # Estimate the thermo for the reactants and products
            item = Reaction(reactants=[m.copy(deep=True) for m in entry.item.reactants], products=[m.copy(deep=True) for m in entry.item.products])
            item.reactants = [Species(molecule=[m]) for m in item.reactants]
            for reactant in item.reactants:
                reactant.generateResonanceIsomers()
                reactant.thermo = thermoDatabase.getThermoData(reactant)
            item.products = [Species(molecule=[m]) for m in item.products]
            for product in item.products:
                product.generateResonanceIsomers()
                product.thermo = thermoDatabase.getThermoData(product)
            # Now that we have the thermo, we can get the reverse k(T)
            item.kinetics = data
            data = item.generateReverseRateCoefficient()
            
            item = Reaction(reactants=entry.item.products, products=entry.item.reactants)
            template = self.getReactionTemplate(item)
            item.degeneracy = self.calculateDegeneracy(item)
            
            new_entry = Entry(
                index = index,
                label = ';'.join([g.label for g in template]),
                item = template,
                data = ArrheniusEP(
                    A = deepcopy(data.A),
                    n = deepcopy(data.n),
                    alpha = 0,
                    E0 = deepcopy(data.Ea),
                    Tmin = deepcopy(data.Tmin),
                    Tmax = deepcopy(data.Tmax),
                ),
                rank = 3,
            )
            new_entry.data.A.value_si /= item.degeneracy
            self.rules.entries[index] = new_entry
            index += 1
    
    def getRootTemplate(self):
        """
        Return the root template for the reaction family. Most of the time this
        is the top-level nodes of the tree (as stored in the 
        :class:`KineticsGroups` object), but there are a few exceptions (e.g.
        R_Recombination).
        """
        if len(self.forwardTemplate.reactants) > len(self.groups.top):
            return self.forwardTemplate.reactants
        else:
            return self.groups.top
    
    def fillKineticsRulesByAveragingUp(self, rootTemplate=None, alreadyDone=None):
        """
        Fill in gaps in the kinetics rate rules by averaging child nodes.
        """
        # If no template is specified, then start at the top-level nodes
        if rootTemplate is None:
            rootTemplate = self.getRootTemplate()
            alreadyDone = {}
        
        rootLabel = ';'.join([g.label for g in rootTemplate])
        
        if rootLabel in alreadyDone:
            return alreadyDone[rootLabel]
        
        if self.hasRateRule(rootTemplate):
            # We already have a rate rule for this exact template
            entry = self.getRateRule(rootTemplate)
            if entry.rank > 0:
                # If the entry has rank of zero, then we have so little faith
                # in it that we'd rather use an averaged value if possible
                # Since this entry does not have a rank of zero, we keep its
                # value
                alreadyDone[rootLabel] = entry.data
                return entry.data
        
        # Recursively descend to the child nodes
        childrenList = [[group] for group in rootTemplate]
        for group in childrenList:
            parent = group.pop(0)
            if len(parent.children) > 0:
                group.extend(parent.children)
            else:
                group.append(parent)
                
        childrenList = getAllCombinations(childrenList)
        kineticsList = []
        for template in childrenList:
            label = ';'.join([g.label for g in template])
            if template == rootTemplate: 
                continue
            
            if label in alreadyDone:
                kinetics = alreadyDone[label]
            else:
                kinetics = self.fillKineticsRulesByAveragingUp(template, alreadyDone)
            
            if kinetics is not None:
                kineticsList.append([kinetics, template])
        
        if len(kineticsList) > 0:
            
            # We found one or more results! Let's average them together
            kinetics = self.__getAverageKinetics([k for k, t in kineticsList])
            kinetics.comment += '(Average of {0})'.format(
                ' + '.join([k.comment if k.comment != '' else ';'.join([g.label for g in t]) for k, t in kineticsList]),
            )
            entry = Entry(
                index = 0,
                label = rootLabel,
                item = rootTemplate,
                data = kinetics,
                rank = 10, # Indicates this is an averaged estimate
            )
            self.rules.entries[entry.label] = entry
            alreadyDone[rootLabel] = entry.data
            return entry.data
            
        alreadyDone[rootLabel] = None
        return None
            
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
        if isinstance(reactantStructures[0], Group):
            reactantStructure = Group()
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
                raise Exception('Unable to change labels from "*" to "*1" and "*2" for reaction family {0}.'.format(label))

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
                        atomLabels['*{0:d}'.format(i)].label = '*{0:d}'.format(4+highest-i)

        if not forward: template = self.reverseTemplate
        else:           template = self.forwardTemplate

        # Split product structure into multiple species if necessary
        productStructures = productStructure.split()
        for product in productStructures:
            product.updateConnectivityValues()

        # Make sure we've made the expected number of products
        if len(template.products) != len(productStructures):
            # We have a different number of products than expected by the template.
            # By definition this means that the template is not a match, so
            # we return None to indicate that we could not generate the product
            # structures
            # We need to think this way in order to distinguish between
            # intermolecular and intramolecular versions of reaction families,
            # which will have very different kinetics
            # Unfortunately this may also squash actual errors with malformed
            # reaction templates
            return None

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

    def __generateProductStructures(self, reactantStructures, maps, forward, **options):
        """
        For a given set of `reactantStructures` and a given set of `maps`,
        generate and return the corresponding product structures. The
        `reactantStructures` parameter should be given in the order the
        reactants are stored in the reaction family template. The `maps`
        parameter is a list of mappings of the top-level tree node of each
        *template* reactant to the corresponding *structure*. This function
        returns the product structures.
        """

        if not forward: template = self.reverseTemplate
        else:           template = self.forwardTemplate

        # Clear any previous atom labeling from all reactant structures
        for struct in reactantStructures: struct.clearLabeledAtoms()

        # Tag atoms with labels
        for m in maps:
            for reactantAtom, templateAtom in m.iteritems():
                reactantAtom.label = templateAtom.label

        # Check that reactant structures are allowed in this family
        # If not, then stop
        for struct in reactantStructures:
            if self.isMoleculeForbidden(struct): raise ForbiddenStructureException()

        # Generate the product structures by applying the forward reaction recipe
        try:
            productStructures = self.applyRecipe(reactantStructures, forward=forward)
            if not productStructures: return None
        except InvalidActionError, e:
            logging.error('Unable to apply reaction recipe!')
            logging.error('Reaction family is {0} in {1} direction'.format(self.label, 'forward' if forward else 'reverse'))
            logging.error('Reactant structures are:')
            for struct in reactantStructures:
                logging.error(struct.toAdjacencyList())
            raise

        # If there are two product structures, place the one containing '*1' first
        if len(productStructures) == 2:
            if not productStructures[0].containsLabeledAtom('*1') and \
                productStructures[1].containsLabeledAtom('*1'):
                productStructures.reverse()

        # Apply the generated species constraints (if given)
        if options:
            maxCarbonAtoms = options.get('maximumCarbonAtoms', 1000000)
            maxHydrogenAtoms = options.get('maximumHydrogenAtoms', 1000000)
            maxOxygenAtoms = options.get('maximumOxygenAtoms', 1000000)
            maxNitrogenAtoms = options.get('maximumNitrogenAtoms', 1000000)
            maxSiliconAtoms = options.get('maximumSiliconAtoms', 1000000)
            maxSulfurAtoms = options.get('maximumSulfurAtoms', 1000000)
            maxHeavyAtoms = options.get('maximumHeavyAtoms', 1000000)
            maxRadicals = options.get('maximumRadicalElectrons', 1000000)
            for struct in productStructures:
                H = struct.getNumAtoms('H')
                if struct.getNumAtoms('C') > maxCarbonAtoms:
                    raise ForbiddenStructureException()
                if H > maxHydrogenAtoms:
                    raise ForbiddenStructureException()
                if struct.getNumAtoms('O') > maxOxygenAtoms:
                    raise ForbiddenStructureException()
                if struct.getNumAtoms('N') > maxNitrogenAtoms:
                    raise ForbiddenStructureException()
                if struct.getNumAtoms('Si') > maxSiliconAtoms:
                    raise ForbiddenStructureException()
                if struct.getNumAtoms('S') > maxSulfurAtoms:
                    raise ForbiddenStructureException()
                if len(struct.atoms) - H > maxHeavyAtoms:
                    raise ForbiddenStructureException()
                if struct.getNumberOfRadicalElectrons() > maxRadicals:
                    raise ForbiddenStructureException()

        # Check that product structures are allowed in this family
        # If not, then stop
        for struct in productStructures:
            if self.isMoleculeForbidden(struct): raise ForbiddenStructureException()

        return productStructures

    def isMoleculeForbidden(self, molecule):
        """
        Return ``True`` if the molecule is forbidden in this family, or
        ``False`` otherwise. 
        """
        from rmgpy.data.rmg import database
        if self.forbidden is not None and self.forbidden.isMoleculeForbidden(molecule):
            return True
        if database.forbiddenStructures.isMoleculeForbidden(molecule):
            return True
        return False

    def __createReaction(self, reactants, products, isForward):
        """
        Create and return a new :class:`Reaction` object containing the
        provided `reactants` and `products` as lists of :class:`Molecule`
        objects.
        """

        # Make sure the products are in fact different than the reactants
        if len(reactants) == len(products) == 1:
            if reactants[0].isIsomorphic(products[0]):
                return None
        elif len(reactants) == len(products) == 2:
            if reactants[0].isIsomorphic(products[0]) and reactants[1].isIsomorphic(products[1]):
                return None
            elif reactants[0].isIsomorphic(products[1]) and reactants[1].isIsomorphic(products[0]):
                return None

        # We need to save the reactant and product structures with atom labels so
        # we can generate the kinetics
        # We make copies so the structures aren't trampled on by later actions
        reactants = [reactant.copy(deep=True) for reactant in reactants]
        #products = [product.copy(deep=True) for product in products]
        for reactant in reactants:
            reactant.updateAtomTypes()
            reactant.updateConnectivityValues()
        for product in products:
            product.updateAtomTypes()
            product.updateConnectivityValues()

        # Create and return reaction object
        return Reaction(reactants=reactants, products=products)

    def __matchReactantToTemplate(self, reactant, templateReactant):
        """
        Return ``True`` if the provided reactant matches the provided
        template reactant and ``False`` if not, along with a complete list of the
        mappings.
        """

        if isinstance(templateReactant, list): templateReactant = templateReactant[0]
        struct = templateReactant.item
        
        if isinstance(struct, LogicNode):
            mappings = []
            for child_structure in struct.getPossibleStructures(self.groups.entries):
                mappings.extend(reactant.findSubgraphIsomorphisms(child_structure))
            return mappings
        elif isinstance(struct, Group):
            return reactant.findSubgraphIsomorphisms(struct)

    def generateReactions(self, reactants, **options):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be either single :class:`Molecule` objects
        or lists of same. Does not estimate the kinetics of these reactions
        at this time. Returns a list of :class:`TemplateReaction` objects
        using :class:`Species` objects for both reactants and products. The
        reactions are constructed such that the forward direction is consistent
        with the template of this reaction family.
        """
        reactionList = []
        
        # Forward direction (the direction in which kinetics is defined)
        reactions = self.__generateReactions(reactants, forward=True, **options)
        for rxn in reactions:
            reaction = TemplateReaction(
                reactants = rxn.reactants[:],
                products = rxn.products[:],
                degeneracy = rxn.degeneracy,
                reversible = rxn.reversible,
                family = self,
            )
            reactionList.append(reaction)
        
        reverseReactions = []
        if self.ownReverse:
            # for each reaction, make its reverse reaction and store in a 'reverse' attribute
            for rxn in reactionList:
                reactions = self.__generateReactions(rxn.products, forward=True, **options)
                reactions = filterReactions(rxn.products, rxn.reactants, reactions)
                assert len(reactions) == 1, "Expecting one matching reverse reaction, not {0}. Forward reaction {1!s} : {1!r}".format(len(reactions), rxn)
                reaction = reactions[0]
                reaction = TemplateReaction(
                    reactants = reaction.reactants[:],
                    products = reaction.products[:],
                    degeneracy = reaction.degeneracy,
                    reversible = reaction.reversible,
                    family = self,
                )
                rxn.reverse = reaction
                reverseReactions.append(reaction)
            
        else: # family is not ownReverse
            # Reverse direction (the direction in which kinetics is not defined)
            reactions = self.__generateReactions(reactants, forward=False, **options)
            for rxn in reactions:
                reaction = TemplateReaction(
                    reactants = rxn.products[:],
                    products = rxn.reactants[:],
                    reversible = rxn.reversible,
                    family = self,
                )
                reaction.degeneracy = self.calculateDegeneracy(reaction)
                reactionList.append(reaction)

        # Determine the reactant-product pairs to use for flux analysis
        # Also store the reaction template (useful so we can easily get the kinetics later)
        for reaction in reactionList:
            reaction.pairs = self.getReactionPairs(reaction)
            reaction.template = self.getReactionTemplate(reaction)
            if hasattr(reaction,'reverse'):
                reaction.reverse.pairs = self.getReactionPairs(reaction.reverse)
                reaction.reverse.template = self.getReactionTemplate(reaction.reverse)
            
        # Return the reactions as containing Species objects, not Molecule objects
        for reaction in reactionList:
            moleculeDict = {}
            for molecule in reaction.reactants:
                moleculeDict[molecule] = Species(molecule=[molecule])
            for molecule in reaction.products:
                moleculeDict[molecule] = Species(molecule=[molecule])
            reaction.reactants = [moleculeDict[molecule] for molecule in reaction.reactants]
            reaction.products = [moleculeDict[molecule] for molecule in reaction.products]
            reaction.pairs = [(moleculeDict[reactant],moleculeDict[product]) for reactant, product in reaction.pairs]

        return reactionList
    
    def calculateDegeneracy(self, reaction):
        """
        For a `reaction` given in the direction in which the kinetics are
        defined, compute the reaction-path degeneracy.
        """
        reactions = self.__generateReactions(reaction.reactants, forward=True)
        products = []
        for product in reaction.products:
            if isinstance(product, Molecule):
                species = Species(molecule=[product])
                species.generateResonanceIsomers()
                products.append(species)
            elif isinstance(product, Species):
               products.append(product.molecule)
        for rxn in reactions:
            # We already know the reactants match, so we only need to evaluate the products
            if len(rxn.products) == len(products) == 1:
                if products[0].isIsomorphic(rxn.products[0]):
                    return rxn.degeneracy
            elif len(rxn.products) == len(products) == 2:
                if products[0].isIsomorphic(rxn.products[0]) and products[1].isIsomorphic(rxn.products[1]):
                    return rxn.degeneracy
                elif products[0].isIsomorphic(rxn.products[1]) and products[1].isIsomorphic(rxn.products[0]):
                    return rxn.degeneracy
        raise Exception('Unable to calculate degeneracy for reaction {0} in reaction family {1}.'.format(reaction, self.label))
    
    def __generateReactions(self, reactants, forward=True, **options):
        """
        Generate a list of all of the possible reactions of this family between
        the list of `reactants`. The number of reactants provided must match
        the number of reactants expected by the template, or this function
        will return an empty list. Each item in the list of reactants should
        be a list of :class:`Molecule` objects, each representing a resonance
        isomer of the species of interest.
        """

        rxnList = []; speciesList = []

        # Wrap each reactant in a list if not already done (this is done to 
        # allow for passing multiple resonance structures for each molecule)
        # This also makes a copy of the reactants list so we don't modify the
        # original
        reactants = [reactant if isinstance(reactant, list) else [reactant] for reactant in reactants]

        sameReactants = len(reactants) == 2 and reactants[0] == reactants[1]
                
        if forward:
            template = self.forwardTemplate
        elif self.reverseTemplate is None:
            return []
        else:
            template = self.reverseTemplate

        # Unimolecular reactants: A --> products
        if len(reactants) == 1 and len(template.reactants) == 1:

            # Iterate over all resonance isomers of the reactant
            for molecule in reactants[0]:

                mappings = self.__matchReactantToTemplate(molecule, template.reactants[0])
                for map in mappings:
                    reactantStructures = [molecule]
                    try:
                        productStructures = self.__generateProductStructures(reactantStructures, [map], forward, **options)
                    except ForbiddenStructureException:
                        pass
                    else:
                        if productStructures is not None:
                            rxn = self.__createReaction(reactantStructures, productStructures, forward)
                            if rxn: rxnList.append(rxn)

        # Bimolecular reactants: A + B --> products
        elif len(reactants) == 2 and len(template.reactants) == 2:

            moleculesA = reactants[0]
            moleculesB = reactants[1]

            # Iterate over all resonance isomers of the reactant
            for moleculeA in moleculesA:
                for moleculeB in moleculesB:

                    # Reactants stored as A + B
                    mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[0])
                    mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[1])

                    # Iterate over each pair of matches (A, B)
                    for mapA in mappingsA:
                        for mapB in mappingsB:
                            reactantStructures = [moleculeA, moleculeB]
                            try:
                                productStructures = self.__generateProductStructures(reactantStructures, [mapA, mapB], forward, **options)
                            except ForbiddenStructureException:
                                pass
                            else:
                                if productStructures is not None:
                                    rxn = self.__createReaction(reactantStructures, productStructures, forward)
                                    if rxn: rxnList.append(rxn)

                    # Only check for swapped reactants if they are different
                    if reactants[0] is not reactants[1]:

                        # Reactants stored as B + A
                        mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[1])
                        mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[0])

                        # Iterate over each pair of matches (A, B)
                        for mapA in mappingsA:
                            for mapB in mappingsB:
                                reactantStructures = [moleculeA, moleculeB]
                                try:
                                    productStructures = self.__generateProductStructures(reactantStructures, [mapA, mapB], forward, **options)
                                except ForbiddenStructureException:
                                    pass
                                else:
                                    if productStructures is not None:
                                        rxn = self.__createReaction(reactantStructures, productStructures, forward)
                                        if rxn: rxnList.append(rxn)

        # Remove duplicates from the reaction list
        index0 = 0
        while index0 < len(rxnList):

            # Generate resonance isomers for products of the current reaction
            products = [product.generateResonanceIsomers() for product in rxnList[index0].products]

            index = index0 + 1
            while index < len(rxnList):
                # We know the reactants are the same, so we only need to compare the products
                match = False
                if len(rxnList[index].products) == len(products) == 1:
                    for product in products[0]:
                        if rxnList[index].products[0].isIsomorphic(product):
                            match = True
                            break
                elif len(rxnList[index].products) == len(products) == 2:
                    for productA in products[0]:
                        for productB in products[1]:
                            if rxnList[index].products[0].isIsomorphic(productA) and rxnList[index].products[1].isIsomorphic(productB):
                                match = True
                                break
                            elif rxnList[index].products[0].isIsomorphic(productB) and rxnList[index].products[1].isIsomorphic(productA):
                                match = True
                                break

                # If we found a match, remove it from the list
                # Also increment the reaction path degeneracy of the remaining reaction
                if match:
                    rxnList.remove(rxnList[index])
                    rxnList[index0].degeneracy += 1
                else:
                    index += 1

            index0 += 1

        # For R_Recombination reactions, the degeneracy is twice what it should
        # be, so divide those by two
        # This is hardcoding of reaction families!
        # For reactions of the form A + A -> products, the degeneracy is twice
        # what it should be, so divide those by two
        if sameReactants or self.label.lower().startswith('r_recombination'):
            for rxn in rxnList:
                assert(rxn.degeneracy % 2 == 0)
                rxn.degeneracy /= 2
                
        # This reaction list has only checked for duplicates within itself, not
        # with the global list of reactions
        return rxnList

    def getReactionPairs(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, return the reactant-product pairs to use when
        performing flux analysis.
        """
        pairs = []; error = False
        if len(reaction.reactants) == 1 or len(reaction.products) == 1:
            # When there is only one reactant (or one product), it is paired 
            # with each of the products (reactants)
            for reactant in reaction.reactants:
                for product in reaction.products:
                    pairs.append([reactant,product])
        elif self.label.lower() == 'h_abstraction':
            # Hardcoding for hydrogen abstraction: pair the reactant containing
            # *1 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*1'):
                if reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
            elif reaction.reactants[1].containsLabeledAtom('*1'):
                if reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
        elif self.label.lower() == 'disproportionation':
            # Hardcoding for disproportionation: pair the reactant containing
            # *1 with the product containing *1
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*1'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
            elif reaction.reactants[1].containsLabeledAtom('*1'):
                if reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
        elif self.label.lower() in ['substitution_o', 'substitutions']:
            # Hardcoding for Substitution_O: pair the reactant containing
            # *2 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*2'):
                if reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
            elif reaction.reactants[1].containsLabeledAtom('*2'):
                if reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
        else:
            error = True
            
        if error:
            raise ReactionPairsError('Unable to determine reaction pairs for {0!s} reaction {1!s}.'.format(self.label, reaction))
        else:
            return pairs
        
    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        return self.groups.getReactionTemplate(reaction)

    def getKineticsForTemplate(self, template, degeneracy=1, method='rate rules'):
        """
        Return an estimate of the kinetics for a reaction with the given
        `template` and reaction-path `degeneracy`. There are two possible methods
        to use: 'group additivity' (new RMG-Py behavior) and 'rate rules' (old
        RMG-Java behavior).
        """
        if method.lower() == 'group additivity':
            return self.estimateKineticsUsingGroupAdditivity(template, degeneracy)
        elif method.lower() == 'rate rules':
            return self.estimateKineticsUsingRateRules(template, degeneracy)
        else:
            raise ValueError('Invalid value "{0}" for method parameter; should be "group additivity" or "rate rules".'.format(method))
        
    def getKineticsFromDepository(self, depository, reaction, template, degeneracy):
        """
        Search the given `depository` in this kinetics family for kinetics
        for the given `reaction`. Returns a list of all of the matching 
        kinetics, the corresponding entries, and ``True`` if the kinetics
        match the forward direction or ``False`` if they match the reverse
        direction.
        """
        kineticsList = []
        if depository.label.endswith('rules'):
            # The depository contains groups
            entries = depository.entries.values()
            for entry in entries:
                entryLabels = entry.label.split(';')
                templateLabels = [group.label for group in template]
                if all([group in entryLabels for group in templateLabels]) and all([group in templateLabels for group in entryLabels]):
                    kineticsList.append([deepcopy(entry.data), entry, True])
            for kinetics, entry, isForward in kineticsList:
                if kinetics is not None:
                    # The rules are defined on a per-site basis, so we need to include the degeneracy manually
                    assert isinstance(kinetics, ArrheniusEP)
                    kinetics.A.value_si *= degeneracy
                    kinetics.comment += "Matched rule {0} {1} in {2}\n".format(entry.index, entry.label, depository.label)
                    kinetics.comment += "Multiplied by reaction path degeneracy {0}".format(degeneracy)
        else:
            # The depository contains real reactions
            entries = depository.entries.values()
            for entry in entries:
                if reaction.isIsomorphic(entry.item):
                    kineticsList.append([deepcopy(entry.data), entry, reaction.isIsomorphic(entry.item, eitherDirection=False)])
            for kinetics, entry, isForward in kineticsList:
                if kinetics is not None:
                    kinetics.comment += "Matched reaction {0} {1} in {2}".format(entry.index, entry.label, depository.label)
        return kineticsList
    
    def getKinetics(self, reaction, template, degeneracy=1, estimator='', returnAllKinetics=True):
        """
        Return the kinetics for the given `reaction` by searching the various
        depositories as well as generating a result using the user-specified `estimator`
        of either 'group additivity' or 'rate rules.'  Unlike
        the regular :meth:`getKinetics()` method, this returns a list of
        results, with each result comprising the kinetics, the source, and
        the entry. If it came from a template estimate, the source and entry
        will both be `None`.
        If returnAllKinetics==False, only the first (best?) matching kinetics is returned.
        """
        kineticsList = []
        
        depositories = self.depositories[:]
        depositories.append(self.rules)
        
        # Check the various depositories for kinetics
        for depository in depositories:
            kineticsList0 = self.getKineticsFromDepository(depository, reaction, template, degeneracy)
            if len(kineticsList0) > 0 and not returnAllKinetics:
                # If we have multiple matching rules but only want one result,
                # choose the one with the lowest rank that occurs first
                if any([x[1].rank == 0 for x in kineticsList0]) and not all([x[1].rank == 0 for x in kineticsList0]):
                    kineticsList0 = [x for x in kineticsList0 if x[1].rank != 0]
                kineticsList0.sort(key=lambda x: (x[1].rank, x[1].index))
                kinetics, entry, isForward = kineticsList0[0]
                return kinetics, depository, entry, isForward
            else:
                for kinetics, entry, isForward in kineticsList0:
                    kineticsList.append([kinetics, depository, entry, isForward])
                    
        # If estimator type of rate rules or group additivity is given, retrieve the kinetics. 
        if estimator:        
            kinetics = self.getKineticsForTemplate(template, degeneracy, method=estimator)
            if kinetics:
                if not returnAllKinetics:
                    return kinetics, None, None, True
                kineticsList.append([kinetics, None, None, True])
        # If no estimation method was given, prioritize rate rule estimation. 
        # If returning all kinetics, add estimations from both rate rules and group additivity.
        else:
            kinetics = self.getKineticsForTemplate(template, degeneracy, method='rate rules')
            if kinetics:
                if not returnAllKinetics:
                    return kinetics, None, None, True
                kineticsList.append([kinetics, 'rate rules', None, True])
            kinetics2 = self.getKineticsForTemplate(template, degeneracy, method='group additivity')
            if kinetics2:
                if not returnAllKinetics:
                    return kinetics, None, None, True
                kineticsList.append([kinetics2, 'group additivity', None, True])
        
        if not returnAllKinetics:
            raise UndeterminableKineticsError(reaction)
        
        return kineticsList
    
    def estimateKineticsUsingGroupAdditivity(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using group additivity.
        """
        # Start with the generic kinetics of the top-level nodes
        kinetics = None
        for entry in self.forwardTemplate.reactants:
            if kinetics is None and entry.data is not None:
                kinetics = entry.data
        if kinetics is None:
            #raise UndeterminableKineticsError('Cannot determine group additivity kinetics estimate for template "{0}".'.format(','.join([e.label for e in template])))
            return None
        # Now add in more specific corrections if possible
        return self.groups.estimateKineticsUsingGroupAdditivity(template, kinetics, degeneracy)
    
    def __getAverageKinetics(self, kineticsList):
        # Although computing via logA is slower, it is necessary because
        # otherwise you could overflow if you are averaging too many values
        logA = 0.0; n = 0.0; E0 = 0.0; alpha = 0.0
        count = len(kineticsList)
        for kinetics in kineticsList:
            logA += math.log10(kinetics.A.value_si)
            n += kinetics.n.value_si
            alpha += kinetics.alpha.value_si
            E0 += kinetics.E0.value_si
        logA /= count
        n /= count
        alpha /= count
        E0 /= count
        Aunits = kineticsList[0].A.units
        if Aunits == 'cm^3/(mol*s)' or 'cm^3/(molecule*s)' or 'm^3/(molecule*s)':
            Aunits = 'm^3/(mol*s)'
        elif Aunits == 'cm^6/(mol^2*s)' or 'cm^6/(molecule^2*s)' or 'm^6/(molecule^2*s)':
            Aunits = 'm^6/(mol^2*s)'
        elif Aunits == 's^-1' or Aunits == 'm^3/(mol*s)' or Aunits == 'm^6/(mol^2*s)':
            pass
        else:
            raise Exception('Invalid units {0} for averaging kinetics.'.format(Aunits))
        averagedKinetics = ArrheniusEP(
            A = (10**logA,Aunits),
            n = n,
            alpha = alpha,
            E0 = (E0*0.001,"kJ/mol"),
        )
        return averagedKinetics
        
        
    def estimateKineticsUsingRateRules(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.
        """
        def getTemplateLabel(template):
            # Get string format of the template in the form "(leaf1,leaf2)"
            return '({0})'.format(','.join([g.label for g in template]))
    
        templateList = [template]
        while len(templateList) > 0:
            
            kineticsList = []
            for t in templateList:
                if self.hasRateRule(t):
                    entry = self.getRateRule(t)
                    kinetics = deepcopy(entry.data)
                    kineticsList.append([kinetics, t])
            
            if len(kineticsList) > 0:                 
                originalLeaves = getTemplateLabel(template)
                                
                if len(kineticsList) == 1:
                    kinetics, t = kineticsList[0]
                    # Check whether the exact rate rule for the original template (most specific
                    # leaves) were found or not.
                    matchedLeaves = getTemplateLabel(t)
                    if matchedLeaves == originalLeaves:
                        kinetics.comment += 'Exact match found' 
                    else:
                    # Using a more general node to estimate original template
                        kinetics.comment += 'Estimated using template ' + matchedLeaves
                else:
                    # We found one or more results! Let's average them together
                    kinetics = self.__getAverageKinetics([k for k, t in kineticsList])
                    kinetics.comment += 'Estimated using average of templates {0}'.format(
                        ' + '.join([getTemplateLabel(t) for k, t in kineticsList]),
                    )
                
                kinetics.comment +=  ' for rate rule ' + originalLeaves
                kinetics.A.value_si *= degeneracy

                return kinetics
            
            else:
                # No results found
                templateList0 = templateList
                templateList = []
                for template0 in templateList0:
                    for index in range(len(template0)):
                        if not template0[index].parent:
                            # We're at the top-level node in this subtreee
                            continue
                        t = template0[:]
                        t[index] = t[index].parent
                        if t not in templateList:
                            templateList.append(t)
                
        # If we're here then we couldn't estimate any kinetics, which is an exception
        raise Exception('Unable to determine kinetics for reaction with template {0}.'.format(template))

################################################################################

def filterReactions(reactants, products, reactionList):
    """
    Remove any reactions from the given `reactionList` whose reactants do
    not involve all the given `reactants` or whose products do not involve 
    all the given `products`. This method checks both forward and reverse
    directions, and only filters out reactions that don't match either.
    """
    
    # Convert from molecules to species and generate resonance isomers.
    reactant_species = []
    for mol in reactants:
        s = Species(molecule=[mol])
        s.generateResonanceIsomers()
        reactant_species.append(s)
    reactants = reactant_species
    product_species = []
    for mol in products:
        s = Species(molecule=[mol])
        s.generateResonanceIsomers()
        product_species.append(s)
    products = product_species
    
    reactions = reactionList[:]
    
    for reaction in reactionList:
        # Forward direction
        reactants0 = [r for r in reaction.reactants]
        for reactant in reactants:
            for reactant0 in reactants0:
                if reactant.isIsomorphic(reactant0):
                    reactants0.remove(reactant0)
                    break
        products0 = [p for p in reaction.products]
        for product in products:
            for product0 in products0:
                if product.isIsomorphic(product0):
                    products0.remove(product0)
                    break
        forward = not (len(reactants0) != 0 or len(products0) != 0)
        # Reverse direction
        reactants0 = [r for r in reaction.products]
        for reactant in reactants:
            for reactant0 in reactants0:
                if reactant.isIsomorphic(reactant0):
                    reactants0.remove(reactant0)
                    break
        products0 = [p for p in reaction.reactants]
        for product in products:
            for product0 in products0:
                if product.isIsomorphic(product0):
                    products0.remove(product0)
                    break
        reverse = not (len(reactants0) != 0 or len(products0) != 0)
        if not forward and not reverse:
            reactions.remove(reaction)
    return reactions
        