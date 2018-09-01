#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains functionality for working with kinetics families.
"""

import os.path
import numpy as np
import logging
import warnings
import codecs
from copy import deepcopy
from collections import OrderedDict

from rmgpy.constraints import failsSpeciesConstraints
from rmgpy.data.base import Database, Entry, LogicNode, LogicOr, ForbiddenStructures,\
                            getAllCombinations
from rmgpy.reaction import Reaction, isomorphic_species_lists
from rmgpy import settings
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Bond, GroupBond, Group, Molecule
from rmgpy.molecule.resonance import generate_aromatic_resonance_structures
from rmgpy.species import Species

from .common import saveEntry, ensure_species, find_degenerate_reactions, generate_molecule_combos,\
                    ensure_independent_atom_ids
from .depository import KineticsDepository
from .groups import KineticsGroups
from .rules import KineticsRules
from rmgpy.exceptions import InvalidActionError, ReactionPairsError, KineticsError,\
                             UndeterminableKineticsError, ForbiddenStructureException,\
                             KekulizationError, ActionError, DatabaseError
import itertools
################################################################################

class TemplateReaction(Reaction):
    """
    A Reaction object generated from a reaction family template. In addition
    to attributes inherited from :class:`Reaction`, this class includes the
    following attributes:

    ============ ========================= =====================================
    Attribute    Type                      Description
    ============ ========================= =====================================
    `family`     ``str``                   The kinetics family that the reaction was created from.
    `estimator`  ``str``                   Whether the kinetics came from rate rules or group additivity.
    `reverse`    :class:`TemplateReaction` The reverse reaction, for families that are their own reverse.
    `is_forward`  ``bool``                 Whether the reaction was generated in the forward direction of the family.
    ============ ========================= =====================================
    """

    def __init__(self,
                index=-1,
                reactants=None,
                products=None,
                specificCollider=None,
                kinetics=None,
                reversible=True,
                transitionState=None,
                duplicate=False,
                degeneracy=1,
                pairs=None,
                family=None,
                template=None,
                estimator=None,
                reverse=None,
                is_forward=None,
                ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          specificCollider=specificCollider,
                          kinetics=kinetics,
                          reversible=reversible,
                          transitionState=transitionState,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs,
                          is_forward=is_forward,
                          )
        self.family = family
        self.template = template
        self.estimator = estimator
        self.reverse = reverse

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TemplateReaction, (self.index,
                                   self.reactants,
                                   self.products,
                                   self.specificCollider,
                                   self.kinetics,
                                   self.reversible,
                                   self.transitionState,
                                   self.duplicate,
                                   self.degeneracy,
                                   self.pairs,
                                   self.family,
                                   self.template,
                                   self.estimator,
                                   self.reverse,
                                   self.is_forward
                                   ))

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        string = 'TemplateReaction('
        if self.index != -1: string += 'index={0:d}, '.format(self.index)
        if self.label != '': string += 'label={0!r}, '.format(self.label)
        if self.reactants is not None: string += 'reactants={0!r}, '.format(self.reactants)
        if self.products is not None: string += 'products={0!r}, '.format(self.products)
        if self.specificCollider is not None: string += 'specificCollider={0!r}, '.format(self.specificCollider)
        if self.kinetics is not None: string += 'kinetics={0!r}, '.format(self.kinetics)
        if not self.reversible: string += 'reversible={0}, '.format(self.reversible)
        if self.transitionState is not None: string += 'transitionState={0!r}, '.format(self.transitionState)
        if self.duplicate: string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1: string += 'degeneracy={0:.1f}, '.format(self.degeneracy)
        if self.pairs is not None: string += 'pairs={0}, '.format(self.pairs)
        if self.family: string += "family='{}', ".format(self.family)
        if self.template: string += "template={}, ".format(self.template)
        if self.comment != '': string += 'comment={0!r}, '.format(self.comment)
        string = string[:-2] + ')'
        return string

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        TemplateReaction this should be a KineticsGroups object.
        """
        return self.family
        
    def copy(self):
        """
        creates a new instance of TemplateReaction
        """
        other = TemplateReaction.__new__(TemplateReaction)
        
        # this was copied from Reaction.copy class
        other.index = self.index
        other.label = self.label
        other.reactants = []
        for reactant in self.reactants:
            other.reactants.append(reactant.copy(deep=True))
        other.products = []
        for product in self.products:
            other.products.append(product.copy(deep=True))
        other.specificCollider = self.specificCollider
        other.degeneracy = self.degeneracy
        other.kinetics = deepcopy(self.kinetics)
        other.reversible = self.reversible
        other.transitionState = deepcopy(self.transitionState)
        other.duplicate = self.duplicate
        other.pairs = deepcopy(self.pairs)
        
        # added for TemplateReaction information
        other.family = self.family
        other.template = self.template
        other.estimator = self.estimator
        other.reverse = self.reverse
        other.is_forward = self.is_forward
        
        return other

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
    GAIN_PAIR     `center`, `pair`              increase the number of lone electron pairs on `center` by `pair`
    LOSE_PAIR     `center`, `pair`              decrease the number of lone electron pairs on `center` by `pair`
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
            elif action[0] == 'LOSE_PAIR':
                other.addAction(['GAIN_PAIR', action[1], action[2]])
            elif action[0] == 'GAIN_PAIR':
                other.addAction(['LOSE_PAIR', action[1], action[2]])
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
        struct.props['validAromatic'] = True

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
                    if bond.isBenzene():
                        struct.props['validAromatic'] = False
                    if doForward:
                        atom1.applyAction(['CHANGE_BOND', label1, info, label2])
                        atom2.applyAction(['CHANGE_BOND', label1, info, label2])
                        bond.applyAction(['CHANGE_BOND', label1, info, label2])
                    else:
                        atom1.applyAction(['CHANGE_BOND', label1, -info, label2])
                        atom2.applyAction(['CHANGE_BOND', label1, -info, label2])
                        bond.applyAction(['CHANGE_BOND', label1, -info, label2])
                elif (action[0] == 'FORM_BOND' and doForward) or (action[0] == 'BREAK_BOND' and not doForward):
                    if struct.hasBond(atom1, atom2):
                        raise InvalidActionError('Attempted to create an existing bond.')
                    bond = GroupBond(atom1, atom2, order=[1]) if pattern else Bond(atom1, atom2, order=1)
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
                        
            elif action[0] in ['LOSE_PAIR', 'GAIN_PAIR']:

                label, change = action[1:]
                change = int(change)

                # Find associated atom
                atom = struct.getLabeledAtom(label)
                if atom is None:
                    raise InvalidActionError('Unable to find atom with label "{0}" while applying reaction recipe.'.format(label))

                # Apply the action
                for i in range(change):
                    if (action[0] == 'GAIN_PAIR' and doForward) or (action[0] == 'LOSE_PAIR' and not doForward):
                        atom.applyAction(['GAIN_PAIR', label, 1])
                    elif (action[0] == 'LOSE_PAIR' and doForward) or (action[0] == 'GAIN_PAIR' and not doForward):
                        atom.applyAction(['LOSE_PAIR', label, 1])

            else:
                raise InvalidActionError('Unknown action "' + action[0] + '" encountered.')

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

class KineticsFamily(Database):
    """
    A class for working with an RMG kinetics family: a set of reactions with 
    similar chemistry, and therefore similar reaction rates. The attributes 
    are:

    =================== =============================== ========================
    Attribute           Type                            Description
    =================== =============================== ========================
    `reverse`           ``string``                      The name of the reverse reaction family
    `reversible`        `Boolean`                       Is family reversible? (True by default)
    `forwardTemplate`   :class:`Reaction`               The forward reaction template
    `forwardRecipe`     :class:`ReactionRecipe`         The steps to take when applying the forward reaction to a set of reactants
    `reverseTemplate`   :class:`Reaction`               The reverse reaction template
    `reverseRecipe`     :class:`ReactionRecipe`         The steps to take when applying the reverse reaction to a set of reactants
    `forbidden`         :class:`ForbiddenStructures`    (Optional) Forbidden product structures in either direction
    `ownReverse`        `Boolean`                       It's its own reverse?
    'boundaryAtoms'     list                            Labels which define the boundaries of end groups in backbone/end families
    `treeDistances`     dict                            The default distance from parent along each tree, if not set default is 1 for every tree
    ------------------- ------------------------------- ------------------------
    `groups`            :class:`KineticsGroups`         The set of kinetics group additivity values
    `rules`             :class:`KineticsRules`          The set of kinetics rate rules from RMG-Java
    `depositories`      ``list``                        A set of additional depositories used to store kinetics data from various sources
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
                 reversible=True,
                 shortDesc='',
                 longDesc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None,
                 boundaryAtoms = None,
                 treeDistances = None
                 ):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
        self.reverse = reverse
        self.reversible = reversible
        self.forwardTemplate = forwardTemplate
        self.forwardRecipe = forwardRecipe
        self.reverseTemplate = reverseTemplate
        self.reverseRecipe = reverseRecipe
        self.forbidden = forbidden
        self.ownReverse = forwardTemplate is not None and reverseTemplate is None
        self.boundaryAtoms = boundaryAtoms
        self.treeDistances = treeDistances
        
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
        warnings.warn("The old kinetics databases are no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
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
            
        self.rules = KineticsRules(label='{0}/rules'.format(self.label))
        self.rules.name = self.rules.label
        try:
            self.rules.loadOld(path, self.groups, numLabels=max(len(self.forwardTemplate.reactants), len(self.groups.top)))
        except Exception:
            logging.error('Error while reading old kinetics family rules from {0!r}.'.format(path))
            raise
        self.depositories = {}

        return self

    def loadOldTemplate(self, path):
        """
        Load an old-style RMG reaction family template from the location `path`.
        """
        warnings.warn("The old kinetics databases are no longer supported and"
                      " may be removed in version 2.3.", DeprecationWarning)
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
        warnings.warn("The old kinetics databases are no longer supported and"
                      " may be removed in version 2.3.", DeprecationWarning)
        if not os.path.exists(path): os.mkdir(path)
        
        self.groups.saveOldDictionary(os.path.join(path, 'dictionary.txt'))
        self.groups.saveOldTree(os.path.join(path, 'tree.txt'))
        # The old kinetics groups use rate rules (not group additivity values),
        # so we can't save the old rateLibrary.txt
        self.saveOldTemplate(os.path.join(path, 'reactionAdjList.txt'))
        # Save forbidden structures if present
        if self.forbidden is not None:
            self.forbidden.saveOld(os.path.join(path, 'forbiddenGroups.txt'))
            
        self.rules.saveOld(path, self)
            
    def saveOldTemplate(self, path):
        """
        Save an old-style RMG reaction family template from the location `path`.
        """
        warnings.warn("The old kinetics databases are no longer supported and"
                      " may be removed in version 2.3.", DeprecationWarning)
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
        
    def distributeTreeDistances(self):
        """
        fills in nodalDistance (the distance between an entry and its parent)
        if not already entered with the value from treeDistances associated
        with the tree the entry comes from
        """
        treeDistances = self.treeDistances
        toplabels = [i.label for i in self.groups.top]
        
        assert len(toplabels) == len(treeDistances), 'treeDistances does not have the same number of entries as there are top nodes in the family'

        for entryName,entry in self.groups.entries.iteritems():
            topentry = entry
            while not (topentry.parent is None): #get the top for the tree entry is in
                topentry = topentry.parent
            if topentry.label in toplabels: #filtering out product nodes
                if entry.nodalDistance is None:
                    entry.nodalDistance = treeDistances[topentry.label]
                
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
        local_context['reverse'] = None
        local_context['reversible'] = None
        local_context['boundaryAtoms'] = None
        local_context['treeDistances'] = None
        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        logging.debug("Loading kinetics family groups from {0}".format(os.path.join(path, 'groups.py')))
        Database.load(self.groups, os.path.join(path, 'groups.py'), local_context, global_context)
        self.name = self.label
        self.boundaryAtoms = local_context.get('boundaryAtoms', None)
        self.treeDistances = local_context.get('treeDistances',None)
        
        # Generate the reverse template if necessary
        self.forwardTemplate.reactants = [self.groups.entries[label] for label in self.forwardTemplate.reactants]
        if self.ownReverse:
            self.forwardTemplate.products = self.forwardTemplate.reactants[:]
            self.reverseTemplate = None
            self.reverseRecipe = None
        else:
            self.reverse = local_context.get('reverse', None)
            self.reversible = True if local_context.get('reversible', None) is None else local_context.get('reversible', None)
            self.forwardTemplate.products = self.generateProductTemplate(self.forwardTemplate.reactants)
            if self.reversible:
                self.reverseTemplate = Reaction(reactants=self.forwardTemplate.products, products=self.forwardTemplate.reactants)
                self.reverseRecipe = self.forwardRecipe.getReverse()
                if self.reverse is None:
                    self.reverse = '{0}_reverse'.format(self.label)
        
        self.groups.numReactants = len(self.forwardTemplate.reactants)
            
        self.rules = KineticsRules(label='{0}/rules'.format(self.label))
        logging.debug("Loading kinetics family rules from {0}".format(os.path.join(path, 'rules.py')))
        self.rules.load(os.path.join(path, 'rules.py'), local_context, global_context)
        
        # load the groups indicated in the entry label
        for label, entries in self.rules.entries.iteritems():
            nodes = label.split(';')
            reactants = [self.groups.entries[node] for node in nodes]
            reaction = Reaction(reactants=reactants, products=[])
            for entry in entries:
                entry.item = reaction
        self.depositories = []
        
        toplabels = [i.label for i in self.groups.top]
        if self.treeDistances is None:
            self.treeDistances = {topentry:1 for topentry in toplabels}

        self.distributeTreeDistances()
            
        if depositoryLabels=='all':
            # Load everything. This option is generally used for working with the database
            # load all the remaining depositories, in order returned by os.walk
            for root, dirs, files in os.walk(path):
                for name in dirs:
                    #if not f.endswith('.py'): continue
                    #name = f.split('.py')[0]
                    #if name not in ['groups', 'rules']:
                    fpath = os.path.join(path, name, 'reactions.py')
                    label = '{0}/{1}'.format(self.label, name)
                    depository = KineticsDepository(label=label)
                    logging.debug("Loading kinetics family depository from {0}".format(fpath))
                    depository.load(fpath, local_context, global_context)
                    self.depositories.append(depository)
                    
            return
                    
        if not depositoryLabels:
            # If depository labels is None or there are no depositories listed, then use the training
            # depository and add them to the RMG rate rules by default:
            depositoryLabels = ['training']
        if depositoryLabels:
            # If there are depository labels, load them in the order specified, but 
            # append the training reactions unless the user specifically declares it not
            # to be included with a '!training' flag
            if '!training' not in depositoryLabels:
                if 'training' not in depositoryLabels:
                    depositoryLabels.append('training')
            
        for name in depositoryLabels :
            if name == '!training':
                continue
            label = '{0}/{1}'.format(self.label, name)
            #f = name+'.py'
            fpath = os.path.join(path, name, 'reactions.py')
            if not os.path.exists(fpath):
                logging.warning("Requested depository {0} does not exist".format(fpath))
                continue
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
            assert action[0] in ['CHANGE_BOND','FORM_BOND','BREAK_BOND','GAIN_RADICAL','LOSE_RADICAL','GAIN_PAIR','LOSE_PAIR']
            self.forwardRecipe.addAction(action)

    def loadForbidden(self, label, group, shortDesc='', longDesc=''):
        """
        Load information about a forbidden structure.
        """
        if not self.forbidden:
            self.forbidden = ForbiddenStructures()
        self.forbidden.loadEntry(label=label, group=group, shortDesc=shortDesc, longDesc=longDesc)

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)
    
    def saveTrainingReactions(self, reactions, reference=None, referenceType='', shortDesc='', longDesc='', rank=3):
        """
        This function takes a list of reactions appends it to the training reactions file.  It ignores the existence of
        duplicate reactions.  
        
        The rank for each new reaction's kinetics is set to a default value of 3 unless the user specifies differently 
        for those reactions.
        
        For each entry, the long description is imported from the kinetics comment. 
        """ 
        from rmgpy import settings

        if not isinstance(reference, list):
            reference = [reference]*len(reactions)
        if not isinstance(referenceType, list):
            referenceType = [referenceType]*len(reactions)
        if not isinstance(shortDesc, list):
            shortDesc = [shortDesc]*len(reactions)
        if not isinstance(longDesc, list):
            longDesc = [longDesc]*len(reactions)
        if not isinstance(rank, list):
            rank = [rank]*len(reactions)

        training_path = os.path.join(settings['database.directory'], 'kinetics', 'families',
                                     self.label, 'training')

        dictionary_path = os.path.join(training_path, 'dictionary.txt')

        # Load the old set of the species of the training reactions
        species_dict = Database().getSpecies(dictionary_path)

        # Add new unique species with labeledAtoms into species_dict
        for rxn in reactions:
            for spec in (rxn.reactants + rxn.products):
                for ex_spec in species_dict.itervalues():
                    if ex_spec.molecule[0].getFormula() != spec.molecule[0].getFormula():
                        continue
                    else:
                        spec_labeled_atoms = spec.molecule[0].getLabeledAtoms()
                        ex_spec_labeled_atoms = ex_spec.molecule[0].getLabeledAtoms()
                        initialMap = {}
                        try:
                            for atomLabel in spec_labeled_atoms:
                                initialMap[spec_labeled_atoms[atomLabel]] = ex_spec_labeled_atoms[atomLabel]
                        except KeyError:
                            # Atom labels did not match, therefore not a match
                            continue
                        if spec.molecule[0].isIsomorphic(ex_spec.molecule[0], initialMap):
                            spec.label = ex_spec.label
                            break
                else:  # No isomorphic existing species found
                    spec_formula = spec.molecule[0].getFormula()
                    if spec_formula not in species_dict:
                        spec.label = spec_formula
                    else:
                        index = 2
                        while (spec_formula + '-{}'.format(index)) in species_dict:
                            index += 1
                        spec.label = spec_formula + '-{}'.format(index)
                    species_dict[spec.label] = spec

        training_file = open(os.path.join(training_path, 'reactions.py'), 'a')

        # get max reaction entry index from the existing training data
        try:
            depository = self.getTrainingDepository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Starting a new one')
            depository = KineticsDepository()
            self.depositories.append(depository)
        
        if depository.entries:
            max_index = max(depository.entries.keys())
        else:
            max_index = 0

        # Add new reactions to training depository
        for i, reaction in enumerate(reactions):    
            index = max_index+i+1
            entry = Entry(
                index = index,
                label = str(reaction),
                item = reaction,
                data = reaction.kinetics,
                reference = reference[i],
                referenceType = referenceType[i],
                shortDesc = unicode(shortDesc[i]),
                longDesc = unicode(longDesc[i]),
                rank = rank[i],
            )

            # Add this entry to the loaded depository so it is immediately usable
            depository.entries[index] = entry
            # Write the entry to the reactions.py file
            self.saveEntry(training_file, entry)

        training_file.close()

        # save species to dictionary
        with open(dictionary_path, 'w') as f:
            for label in species_dict.keys():
                f.write(species_dict[label].molecule[0].toAdjacencyList(label=label, removeH=False))
                f.write('\n')

    def save(self, path):
        """
        Save the current database to the file at location `path` on disk. 
        """
        self.saveGroups(os.path.join(path, 'groups.py'))
        self.rules.save(os.path.join(path, 'rules.py'))
        for depository in self.depositories:
            self.saveDepository(depository, os.path.join(path, '{0}'.format(depository.label[len(self.label)+1:])))
    
    def saveDepository(self, depository, path):
        """
        Save the given kinetics family `depository` to the location `path` on
        disk.
        """
        depository.saveDictionary(os.path.join(path,'dictionary.txt'))
        depository.save(os.path.join(path,'reactions.py'))
        
    def saveGroups(self, path):
        """
        Save the current database to the file at location `path` on disk. 
        """
        entries = self.groups.getEntriesToSave()
                
        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}/groups"\n'.format(self.name))
        f.write('shortDesc = u"{0}"\n'.format(self.groups.shortDesc))
        f.write('longDesc = u"""\n')
        f.write(self.groups.longDesc)
        f.write('\n"""\n\n')

        # Write the template
        f.write('template(reactants=[{0}], products=[{1}], ownReverse={2})\n\n'.format(
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.reactants]),
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.products]),
            self.ownReverse))

        # Write reverse name
        if not self.ownReverse:
            if self.reverse is not None:
                f.write('reverse = "{0}"\n'.format(self.reverse))
            else:
                f.write('reverse = None\n')
        
        f.write('reversible = {0}\n\n'.format(self.reversible))

        # Write the recipe
        f.write('recipe(actions=[\n')
        for action in self.forwardRecipe.actions:
            f.write('    {0!r},\n'.format(action))
        f.write('])\n\n')

        if self.boundaryAtoms:
            f.write('boundaryAtoms = ["{0}", "{1}"]'.format(self.boundaryAtoms[0], self.boundaryAtoms[1]))
            f.write('\n\n')

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
        logging.log(1, "Generating template for products.")
        for reactant in reactants0:
            if isinstance(reactant, list):  reactants = [reactant[0]]
            else:                           reactants = [reactant]

            logging.log(1, "Reactants: {0}".format(reactants))
            for s in reactants:
                logging.log(1, "Reactant {0}".format(s))
                struct = s.item
                if isinstance(struct, LogicNode):
                    all_structures = struct.getPossibleStructures(self.groups.entries)
                    logging.log(1, 'Expanding logic node {0} to {1}'.format(s, all_structures))
                    reactantStructures.append(all_structures)
                    for p in all_structures:
                        logging.log(1, p.toAdjacencyList() )
                else:
                    reactantStructures.append([struct])
                    logging.log(1, struct.toAdjacencyList() )

        # Second, get all possible combinations of reactant structures
        reactantStructures = getAllCombinations(reactantStructures)
        
        # Third, generate all possible product structures by applying the
        # recipe to each combination of reactant structures
        # Note that bimolecular products are split by labeled atoms
        productStructures = []
        for reactantStructure in reactantStructures:
            productStructure = self.applyRecipe(reactantStructure, forward=True, unique=False)
            if productStructure:
                productStructures.append(productStructure)

        # Fourth, remove duplicates from the lists
        productStructureList = [[] for i in range(len(productStructures[0]))]
        for productStructure in productStructures:
            for i, struct in enumerate(productStructure):
                for s in productStructureList[i]:
                    try:
                        if s.isIdentical(struct): break
                    except KeyError:
                        logging.error(struct.toAdjacencyList())
                        logging.error(s.toAdjacencyList())
                        raise
                else:
                    productStructureList[i].append(struct)
                    
        logging.log(1, "Unique generated product structures:")
        logging.log(1, "\n".join([p[0].toAdjacencyList() for p in productStructures]))
        
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
                children = []
                counter = 0
                for product in products:
                    entry = Entry(
                        label = '{0}{1:d}'.format(label,counter+1),
                        item = product,
                    )                
                    children.append(entry)
                    self.groups.entries[entry.label] = entry
                    counter += 1
                
                # Enter the parent of the groups as a logicOr of all the products
                entry = Entry(
                    label = label,
                    item = LogicOr([child.label for child in children],invert=False),
                    children = children,
                )
                self.groups.entries[entry.label] = entry
                # Make this entry the parent of all its children
                for child in children:
                    child.parent = entry
                counter += 1
                productSet.append(entry)

        return productSet

    def hasRateRule(self, template):
        """
        Return ``True`` if a rate rule with the given `template` currently 
        exists, or ``False`` otherwise.
        """
        return self.rules.hasRule(template)

    def getRateRule(self, template):
        """
        Return the rate rule with the given `template`. Raises a 
        :class:`ValueError` if no corresponding entry exists.
        """
        entry = self.rules.getRule(template)
        if entry is None:
            raise ValueError('No entry for template {0}.'.format(template))
        return entry

    def addKineticsRulesFromTrainingSet(self, thermoDatabase=None):
        """
        For each reaction involving real reactants and products in the training
        set, add a rate rule for that reaction.
        """
        try:
            depository = self.getTrainingDepository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Must be because you turned off the training depository.')
            return
        
        tentries = depository.entries
        
        index = max([e.index for e in self.rules.getEntries()] or [0]) + 1
        
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
            
            tentries[entry.index].item.is_forward = True
            
            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            
            new_entry = Entry(
                index = index,
                label = ';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template],
                                                   products=[]),
                data = data.toArrheniusEP(),
                rank = entry.rank,
                reference=entry.reference,
                shortDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.shortDesc,
                longDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.longDesc,
            )
            new_entry.data.comment = "From training reaction {1} used for {0}".format(';'.join([g.label for g in template]), entry.index)

            new_entry.data.A.value_si /= entry.item.degeneracy
            try:
                self.rules.entries[new_entry.label].append(new_entry)
            except KeyError:
                self.rules.entries[new_entry.label] = [new_entry]
            index += 1
        
        # Process the entries that are stored in the reverse direction of the
        # family definition
        for entry in reverse_entries:
            
            tentries[entry.index].item.is_forward = False
            
            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            # Estimate the thermo for the reactants and products
            # trainingSet=True used later to does not allow species to match a liquid phase library and get corrected thermo which will affect reverse rate calculation
            item = Reaction(reactants=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label) for m in entry.item.reactants],
                             products=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label) for m in entry.item.products])
            for reactant in item.reactants:
                reactant.generate_resonance_structures()
                reactant.thermo = thermoDatabase.getThermoData(reactant, trainingSet=True)
            for product in item.products:
                product.generate_resonance_structures()
                product.thermo = thermoDatabase.getThermoData(product,trainingSet=True)
            # Now that we have the thermo, we can get the reverse k(T)
            item.kinetics = data
            data = item.generateReverseRateCoefficient()
            
            item = TemplateReaction(reactants=[m.molecule[0].copy(deep=True) for m in entry.item.products], 
                                               products=[m.molecule[0].copy(deep=True) for m in entry.item.reactants])
            template = self.getReactionTemplate(item)

            item.template = self.getReactionTemplateLabels(item)
            new_degeneracy = self.calculateDegeneracy(item)

            new_entry = Entry(
                index = index,
                label = ';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template],
                                                   products=[]),
                data = data.toArrheniusEP(),
                rank = entry.rank,
                reference=entry.reference,
                shortDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.shortDesc,
                longDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.longDesc,
            )
            new_entry.data.comment = "From training reaction {1} used for {0}".format(';'.join([g.label for g in template]), entry.index)

            new_entry.data.A.value_si /= new_degeneracy
            try:
                self.rules.entries[new_entry.label].append(new_entry)
            except KeyError:
                self.rules.entries[new_entry.label] = [new_entry]
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
    
    def fillKineticsRulesByAveragingUp(self, verbose=False):
        """
        Fill in gaps in the kinetics rate rules by averaging child nodes
        recursively starting from the top level root template.
        """
        
        self.rules.fillRulesByAveragingUp(self.getRootTemplate(), {}, verbose)

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

        if forward:
            # Hardcoding of reaction family for radical recombination (colligation)
            # because the two reactants are identical, they have the same tags
            # In this case, we must change the labels from '*' and '*' to '*1' and
            # '*2'
            if label == 'r_recombination':
                identicalCenterCounter = 0
                for atom in reactantStructure.atoms:
                    if atom.label == '*':
                        identicalCenterCounter += 1
                        atom.label = '*' + str(identicalCenterCounter)
                if identicalCenterCounter != 2:
                    raise KineticsError(
                        'Trying to apply recipe for reaction family {}: Only one occurrence of "*" found.'.format(label)
                    )
            # Hardcoding of reaction family for peroxyl disproportionation
            # '*1' and '*2' have to be changed to '*3' and '*4' for the second reactant
            elif label == 'peroxyl_disproportionation':
                identicalCenterCounter1 = identicalCenterCounter2 = 0
                for atom in reactantStructure.atoms:
                    if atom.label == '*1':
                        identicalCenterCounter1 += 1
                        if identicalCenterCounter1 > 1:
                            atom.label = '*3'
                    elif atom.label == '*2':
                        identicalCenterCounter2 += 1
                        if identicalCenterCounter2 > 1:
                            atom.label = '*4'
                msg = 'Trying to apply recipe for reaction family {}:'.format(label)
                error = False
                if identicalCenterCounter1 != 2:
                    msg += ' Only one occurrence of "*1" found.'
                    error = True
                if identicalCenterCounter2 != 2:
                    msg += ' Only one occurrence of "*2" found.'
                    error = True
                if error:
                    raise KineticsError(msg)
            # Hardcoding of reaction family for bimolecular hydroperoxide decomposition
            # '*2' has to be changed to '*4' for the second reactant and '*1' has to be
            # changed to '*6'. '*3' has to be changed to '*5' for the first reactant.
            # '*5' and '*6' do no participate in the reaction but are required for
            # relabeling in the reverse direction.
            elif label == 'bimolec_hydroperoxide_decomposition':
                identicalCenterCounter1 = identicalCenterCounter2 = identicalCenterCounter3 = 0
                for atom in reactantStructure.atoms:
                    if atom.label == '*1':
                        identicalCenterCounter1 += 1
                        if identicalCenterCounter1 > 1:
                            atom.label = '*6'
                    elif atom.label == '*2':
                        identicalCenterCounter2 += 1
                        if identicalCenterCounter2 > 1:
                            atom.label = '*4'
                    elif atom.label == '*3':
                        identicalCenterCounter3 += 1
                        if identicalCenterCounter3 == 1:
                            atom.label = '*5'
                msg = 'Trying to apply recipe for reaction family {}:'.format(label)
                error = False
                if identicalCenterCounter1 != 2:
                    msg += ' Only one occurrence of "*1" found.'
                    error = True
                if identicalCenterCounter2 != 2:
                    msg += ' Only one occurrence of "*2" found.'
                    error = True
                if identicalCenterCounter3 != 2:
                    msg += ' Only one occurrence of "*3" found.'
                    error = True
                if error:
                    raise KineticsError(msg)

            # Generate the product structure by applying the recipe
            self.forwardRecipe.applyForward(reactantStructure, unique)
        else:
            self.reverseRecipe.applyForward(reactantStructure, unique)

        if not reactantStructure.props['validAromatic']:
            if isinstance(reactantStructure, Molecule):
                # For molecules, kekulize the product to redistribute bonds appropriately
                reactantStructure.kekulize()
            else:
                # For groups, we ignore the product template for a purely aromatic group
                # If there is an analagous aliphatic group in the family, then the product template will be identical
                # There should NOT be any families that consist solely of aromatic reactant templates
                return []
        productStructure = reactantStructure

        if not forward:
            # Hardcoding of reaction family for reverse of radical recombination
            # (Unimolecular homolysis)
            # Because the two products are identical, they should the same tags
            # In this case, we must change the labels from '*1' and '*2' to '*' and
            # '*'
            if label == 'r_recombination':
                for atom in productStructure.atoms:
                    if atom.label == '*1' or atom.label == '*2':
                        atom.label = '*'
            # Hardcoding of reaction family for reverse of peroxyl disproportionation
            # Labels '*3' and '*4' have to be changed back to '*1' and '*2'
            elif label == 'peroxyl_disproportionation':
                for atom in productStructure.atoms:
                    if atom.label == '*3':
                        atom.label = '*1'
                    elif atom.label == '*4':
                        atom.label = '*2'
            # Hardcoding of reaction family for bimolecular hydroperoxide decomposition
            # '*5' has to be changed back to '*3', '*6' has to be changed to '*1', and
            # '*4' has to be changed to '*2'
            elif label == 'bimolec_hydroperoxide_decomposition':
                for atom in productStructure.atoms:
                    if atom.label == '*5':
                        atom.label = '*3'
                    elif atom.label == '*6':
                        atom.label = '*1'
                    elif atom.label == '*4':
                        atom.label = '*2'

        # If reaction family is its own reverse, relabel atoms
        # This allows comparison of the product species to forbidden
        #  structures which are labeled as reactants.
        # Unfortunately, this means that reaction family info is
        #  hardcoded, so this must be updated if the database changes.
        if not self.reverseTemplate:
            # Get atom labels for products
            atomLabels = {}
            for atom in productStructure.atoms:
                if atom.label != '':
                    atomLabels[atom.label] = atom

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
                highest = len(atomLabels)
                if highest > 4:
                    # swap *4 with *5
                    atomLabels['*4'].label = '*5'
                    atomLabels['*5'].label = '*4'
                if highest > 6:
                    # swap *6 with the highest, etc.
                    for i in range(6, highest+1):
                        atomLabels['*{0:d}'.format(i)].label = '*{0:d}'.format(6+highest-i)
                        
            elif label == 'intra_ene_reaction':
                # Labels for nodes are swapped
                atomLabels['*1'].label = '*2'
                atomLabels['*2'].label = '*1'
                atomLabels['*3'].label = '*5'
                atomLabels['*5'].label = '*3'

            elif label == '6_membered_central_c-c_shift':
                # Labels for nodes are swapped
                atomLabels['*1'].label = '*3'
                atomLabels['*3'].label = '*1'
                atomLabels['*4'].label = '*6'
                atomLabels['*6'].label = '*4'

            elif label == '1,2_shiftc':
                # Labels for nodes are swapped
                atomLabels['*2'].label = '*3'
                atomLabels['*3'].label = '*2'

            elif label == 'intra_r_add_exo_scission':
                # Labels for nodes are swapped
                atomLabels['*1'].label = '*3'
                atomLabels['*3'].label = '*1'

            elif label == 'intra_substitutions_isomerization':
                # Swap *2 and *3
                atomLabels['*2'].label = '*3'
                atomLabels['*3'].label = '*2'

        if not forward:
            template = self.reverseTemplate
        else:
            template = self.forwardTemplate

        # Split product structure into multiple species if necessary
        productStructures = productStructure.split()

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

        # Make sure we don't create a different net charge between reactants and products
        reactant_net_charge = product_net_charge = 0
        for struc in reactantStructures:
            struc.update()
            reactant_net_charge += struc.getNetCharge()
        for struct in productStructures:
            # If product structures are Molecule objects, update their atom types
            # If product structures are Group objects and the reaction is in certain families
            # (families with charged substances), the charge of structures will be updated
            if isinstance(struct, Molecule):
                struct.update()
            elif isinstance(struct, Group):
                struct.resetRingMembership()
                if label in ['1,2_insertion_co', 'r_addition_com', 'co_disproportionation',
                             'intra_no2_ono_conversion', 'lone_electron_pair_bond']:
                    struct.update_charge()
            else:
                raise TypeError('Expecting Molecule or Group object, not {0}'.format(struct.__class__.__name__))
            product_net_charge += struc.getNetCharge()
        if reactant_net_charge != product_net_charge:
            logging.debug('The net charge of the reactants {0} differs from the net charge of the products {1} in'
                          ' reaction family {2}. Not generating this reaction.'.format(
                           reactant_net_charge,product_net_charge,self.label))
            return None
        # The following check should be removed once RMG can process charged species
        # This is applied only for :class:Molecule (not for :class:Group which is allowed to have a nonzero net charge)
        if any([structure.getNetCharge() for structure in reactantStructures + productStructures])\
                and isinstance(struc, Molecule):
            logging.debug('A net charged species was formed when reacting {0} to form {1} in'
                          ' reaction family {2}. Not generating this reaction.'.format(
                           reactant_net_charge,product_net_charge,self.label))
            return None

        # If there are two product structures, place the one containing '*1' first
        if len(productStructures) == 2:
            if not productStructures[0].containsLabeledAtom('*1') and\
                    productStructures[1].containsLabeledAtom('*1'):
                productStructures.reverse()
        # If there are three product structures, sort them based on the lowest number label in each structure
        elif len(productStructures) == 3:
            lowest_labels = []
            for struct in productStructures:
                # Extract digits from labels and convert others (e.g., "*") to empty strings
                labels = [''.join(c for c in label if c.isdigit()) for label in struct.getLabeledAtoms().keys()]
                # Convert digits to integers and remove empty strings
                labels = [int(label) for label in labels if label]
                lowest_labels.append(min(labels))
            productStructures = [s for _, s in sorted(zip(lowest_labels, productStructures))]
            
        # Return the product structures
        return productStructures

    def __generateProductStructures(self, reactantStructures, maps, forward):
        """
        For a given set of `reactantStructures` and a given set of `maps`,
        generate and return the corresponding product structures. The
        `reactantStructures` parameter should be given in the order the
        reactants are stored in the reaction family template. The `maps`
        parameter is a list of mappings of the top-level tree node of each
        *template* reactant to the corresponding *structure*. This function
        returns a list of the product structures.
        """
        
        # Clear any previous atom labeling from all reactant structures
        for struct in reactantStructures: struct.clearLabeledAtoms()

        # Tag atoms with labels
        for m in maps:
            for reactantAtom, templateAtom in m.iteritems():
                reactantAtom.label = templateAtom.label

        # Check that reactant structures are allowed in this family
        # If not, then stop
        for struct in reactantStructures:
            if self.isMoleculeForbidden(struct):
                raise ForbiddenStructureException()

        # Generate the product structures by applying the forward reaction recipe
        try:
            productStructures = self.applyRecipe(reactantStructures, forward=forward)
            if not productStructures: return None
        except (InvalidActionError, KekulizationError):
            # If unable to apply the reaction recipe, then return no product structures
            return None
        except ActionError:
            logging.error(
                'Could not generate product structures for reaction family {0} in {1} direction'.format(
                    self.label, 'forward' if forward else 'reverse'))
            logging.info('Reactant structures:')
            for struct in reactantStructures:
                logging.info('{0}\n{1}\n'.format(struct, struct.toAdjacencyList()))
            raise

        # Apply the generated species constraints (if given)
        for struct in productStructures:
            if self.isMoleculeForbidden(struct):
                raise ForbiddenStructureException() 
            if failsSpeciesConstraints(struct):
                raise ForbiddenStructureException() 
                
        return productStructures

    def isMoleculeForbidden(self, molecule):
        """
        Return ``True`` if the molecule is forbidden in this family, or
        ``False`` otherwise. 
        """

        # check family-specific forbidden structures 
        if self.forbidden is not None and self.forbidden.isMoleculeForbidden(molecule):
            return True


        return False

    def __createReaction(self, reactants, products, is_forward):
        """
        Create and return a new :class:`Reaction` object containing the
        provided `reactants` and `products` as lists of :class:`Molecule`
        objects.
        """

        # Make sure the products are in fact different than the reactants
        if isomorphic_species_lists(reactants, products):
            return None

        # Create and return template reaction object
        reaction = TemplateReaction(
            reactants = reactants if is_forward else products,
            products = products if is_forward else reactants,
            degeneracy = 1,
            reversible = self.reversible,
            family = self.label,
            is_forward = is_forward,
        )
        
        # Store the labeled atoms so we can recover them later
        # (e.g. for generating reaction pairs and templates)
        labeledAtoms = []
        for reactant in reaction.reactants:
            for label, atom in reactant.getLabeledAtoms().items():
                labeledAtoms.append((label, atom))
        reaction.labeledAtoms = labeledAtoms
        
        return reaction

    def __matchReactantToTemplate(self, reactant, templateReactant):
        """
        Return a complete list of the mappings if the provided reactant 
        matches the provided template reactant, or an empty list if not.
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
        else:
            raise NotImplementedError("Not expecting template of type {}".format(type(struct)))

    def generateReactions(self, reactants, products=None, prod_resonance=True):
        """
        Generate all reactions between the provided list of one, two, or three
        `reactants`, which should be either single :class:`Molecule` objects
        or lists of same. Does not estimate the kinetics of these reactions
        at this time. Returns a list of :class:`TemplateReaction` objects
        using :class:`Molecule` objects for both reactants and products
        The reactions are constructed such that the forward direction is
        consistent with the template of this reaction family.

        Args:
            reactants (list):                List of Molecules to react.
            products (list, optional):       List of Molecules or Species of desired product structures.
            prod_resonance (bool, optional): Flag to generate resonance structures for product checking.
                Defaults to True, resonance structures are compared.

        Returns:
            List of all reactions containing Molecule objects with the
            specified reactants and products within this family.
            Degenerate reactions are returned as separate reactions.
        """
        reactionList = []

        # Forward direction (the direction in which kinetics is defined)
        reactionList.extend(
            self.__generateReactions(reactants, products=products, forward=True, prod_resonance=prod_resonance))

        if not self.ownReverse and self.reversible:
            # Reverse direction (the direction in which kinetics is not defined)
            reactionList.extend(
                self.__generateReactions(reactants, products=products, forward=False, prod_resonance=prod_resonance))

        return reactionList

    def addReverseAttribute(self, rxn, react_non_reactive=True):
        """
        For rxn (with species' objects) from families with ownReverse, this method adds a `reverse`
        attribute that contains the reverse reaction information (like degeneracy)

        Returns `True` if successful and `False` if the reverse reaction is forbidden.
        Will raise a `KineticsError` if unsuccessful for other reasons.
        """
        if self.ownReverse and all([spc.has_reactive_molecule() for spc in rxn.products]):
            # Check if the reactants are the same
            sameReactants = 0
            if len(rxn.products) == 2 and rxn.products[0].isIsomorphic(rxn.products[1]):
                sameReactants = 2
            elif len(rxn.products) == 3:
                same_01 = rxn.products[0].isIsomorphic(rxn.products[1])
                same_02 = rxn.products[0].isIsomorphic(rxn.products[2])
                if same_01 and same_02:
                    sameReactants = 3
                elif same_01 or same_02:
                    sameReactants = 2
                elif rxn.products[1].isIsomorphic(rxn.products[2]):
                    sameReactants = 2

            reactionList = self.__generateReactions([spc.molecule for spc in rxn.products],
                                                    products=rxn.reactants, forward=True,
                                                    react_non_reactive=react_non_reactive)
            reactions = find_degenerate_reactions(reactionList, sameReactants, kinetics_family=self)
            if len(reactions) == 0:
                logging.error("Expecting one matching reverse reaction, not zero in reaction family {0} for forward reaction {1}.\n".format(self.label, str(rxn)))
                logging.error("There is likely a bug in the RMG-database kinetics reaction family involving a missing group, missing atomlabels, forbidden groups, etc.")
                for reactant in rxn.reactants:
                    logging.info("Reactant")
                    logging.info(reactant.toAdjacencyList())
                for product in rxn.products:
                    logging.info("Product")
                    logging.info(product.toAdjacencyList())
                logging.error("Debugging why no reaction was found...")
                logging.error("Checking whether the family's forbidden species have affected reaction generation...")
                # Set family's forbidden structures to empty for now to see if reaction gets generated...
                # Note that it is not necessary to check global forbidden structures, because this reaction would not have
                # been formed in the first place.
                tempObject = self.forbidden
                self.forbidden = ForbiddenStructures()  # Initialize with empty one
                try:
                    reactionList = self.__generateReactions([spc.molecule for spc in rxn.products],
                                                            products=rxn.reactants, forward=True,
                                                            react_non_reactive=react_non_reactive)
                    reactions = find_degenerate_reactions(reactionList, sameReactants, kinetics_family=self)
                finally:
                    self.forbidden = tempObject
                if len(reactions) == 1 or (len(reactions) > 1 and all([reactions[0].isIsomorphic(other, checkTemplateRxnProducts=True) for other in reactions])):
                    logging.error("Error was fixed, the product is a forbidden structure when used as a reactant in the reverse direction.")
                    # This reaction should be forbidden in the forward direction as well
                    return False
                else:
                    logging.error("Still experiencing error: Expecting one matching reverse reaction, not {0} in reaction family {1} for forward reaction {2}.\n".format(len(reactions), self.label, str(rxn)))
                    raise KineticsError("Did not find reverse reaction in reaction family {0} for reaction {1}.".format(self.label, str(rxn)))
            elif len(reactions) > 1 and not all([reactions[0].isIsomorphic(other, checkTemplateRxnProducts=True) for other in reactions]):
                logging.error("Expecting one matching reverse reaction. Recieved {0} reactions with multiple non-isomorphic ones in reaction family {1} for forward reaction {2}.\n".format(len(reactions), self.label, str(rxn)))
                logging.info("Found the following reverse reactions")
                for rxn0 in reactions:
                    logging.info(str(rxn0))
                    for reactant in rxn0.reactants:
                        logging.info("Reactant")
                        logging.info(reactant.toAdjacencyList())
                    for product in rxn0.products:
                        logging.info("Product")
                        logging.info(product.toAdjacencyList())
                raise KineticsError("Found multiple reverse reactions in reaction family {0} for reaction {1}, likely due to inconsistent resonance structure generation".format(self.label, str(rxn)))
            else:
                rxn.reverse = reactions[0]
                return True

    def calculateDegeneracy(self, reaction):
        """
        For a `reaction`  with `Molecule` or `Species` objects given in the direction in which
        the kinetics are defined, compute the reaction-path degeneracy.

        This method by default adjusts for double counting of identical reactants. 
        This should only be adjusted once per reaction. To not adjust for 
        identical reactants (since you will be reducing them later in the algorithm), add
        `ignoreSameReactants= True` to this method.
        """
        # Check if the reactants are the same
        # If they refer to the same memory address, then make a deep copy so
        # they can be manipulated independently
        reactants = reaction.reactants
        same_reactants = 0
        if len(reactants) == 2:
            if reactants[0] is reactants[1]:
                reactants[1] = reactants[1].copy(deep=True)
                same_reactants = 2
            elif reactants[0].isIsomorphic(reactants[1]):
                same_reactants = 2
        elif len(reactants) == 3:
            same_01 = reactants[0] is reactants[1]
            same_02 = reactants[0] is reactants[2]
            if same_01 and same_02:
                same_reactants = 3
                reactants[1] = reactants[1].copy(deep=True)
                reactants[2] = reactants[2].copy(deep=True)
            elif same_01:
                same_reactants = 2
                reactants[1] = reactants[1].copy(deep=True)
            elif same_02:
                same_reactants = 2
                reactants[2] = reactants[2].copy(deep=True)
            elif reactants[1] is reactants[2]:
                same_reactants = 2
                reactants[2] = reactants[2].copy(deep=True)
            else:
                same_01 = reactants[0].isIsomorphic(reactants[1])
                same_02 = reactants[0].isIsomorphic(reactants[2])
                if same_01 and same_02:
                    same_reactants = 3
                elif same_01 or same_02:
                    same_reactants = 2
                elif reactants[1].isIsomorphic(reactants[2]):
                    same_reactants = 2

        # Label reactant atoms for proper degeneracy calculation
        ensure_independent_atom_ids(reactants, resonance=True)
        molecule_combos = generate_molecule_combos(reactants)

        reactions = []
        for combo in molecule_combos:
            reactions.extend(self.__generateReactions(combo, products=reaction.products, forward=True,
                                                      react_non_reactive=True))

        # remove degenerate reactions
        reactions = find_degenerate_reactions(reactions, same_reactants, template=reaction.template, kinetics_family=self)

        # log issues
        if len(reactions) != 1:
            for reactant in reaction.reactants:
                logging.error("Reactant: {0!r}".format(reactant))
            for product in reaction.products:
                logging.error("Product: {0!r}".format(product))
            raise KineticsError(('Unable to calculate degeneracy for reaction {0} '
                                 'in reaction family {1}. Expected 1 reaction '
                                 'but generated {2}').format(reaction, self.label, len(reactions)))
        return reactions[0].degeneracy
        
    def __generateReactions(self, reactants, products=None, forward=True, prod_resonance=True,
                            react_non_reactive=False):
        """
        Generate a list of all the possible reactions of this family between
        the list of `reactants`. The number of reactants provided must match
        the number of reactants expected by the template, or this function
        will return an empty list. Each item in the list of reactants should
        be a list of :class:`Molecule` objects, each representing a resonance
        structure of the species of interest.
        
        This method returns all reactions, and degenerate reactions can then be
        found using `rmgpy.data.kinetics.common.find_degenerate_reactions`.

        Args:
            reactants:          List of Molecules to react
            products:           List of Molecules or Species of desired product structures (optional)
            forward:            Flag to indicate whether the forward or reverse template should be applied (optional)
                                Default is True, forward template is used
            prod_resonance:     Flag to generate resonance structures for product checking (optional)
                                Default is True, resonance structures are compared
            react_non_reactive: Flag to generate reactions between unreactive molecules (optional)
                                Default is False, reactions involving unreactive molecules are not generated

        Returns:
            List of all reactions containing Molecule objects with the
                specified reactants and products within this family.
            Degenerate reactions are returned as separate reactions.
        """

        rxnList = []; speciesList = []

        # Wrap each reactant in a list if not already done (this is done to 
        # allow for passing multiple resonance structures for each molecule)
        # This also makes a copy of the reactants list so we don't modify the
        # original
        reactants = [reactant if isinstance(reactant, list) else [reactant] for reactant in reactants]

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
                if molecule.reactive or react_non_reactive:  # don't react non representative resonance isomers unless
                    # explicitly desired (e.g., when called from calculateDegeneracy)
                    mappings = self.__matchReactantToTemplate(molecule, template.reactants[0])
                    for map in mappings:
                        reactantStructures = [molecule]
                        try:
                            productStructures = self.__generateProductStructures(reactantStructures, [map], forward)
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
                    if (moleculeA.reactive and moleculeB.reactive) or react_non_reactive:

                        # Reactants stored as A + B
                        mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[0])
                        mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[1])

                        # Iterate over each pair of matches (A, B)
                        for mapA in mappingsA:
                            for mapB in mappingsB:
                                # Reverse the order of reactants in case we have a family with only one reactant tree
                                # that can produce different products depending on the order of reactants
                                reactantStructures = [moleculeB, moleculeA]
                                try:
                                    productStructures = self.__generateProductStructures(reactantStructures, [mapB, mapA], forward)
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
                                        productStructures = self.__generateProductStructures(reactantStructures, [mapA, mapB], forward)
                                    except ForbiddenStructureException:
                                        pass
                                    else:
                                        if productStructures is not None:
                                            rxn = self.__createReaction(reactantStructures, productStructures, forward)
                                            if rxn: rxnList.append(rxn)
        
        # Trimolecular reactants: A + B + C --> products
        elif len(reactants) == 3 and len(template.reactants) == 3:

            moleculesA = reactants[0]
            moleculesB = reactants[1]
            moleculesC = reactants[2]

            # Iterate over all resonance isomers of the reactants
            for moleculeA in moleculesA:
                for moleculeB in moleculesB:
                    for moleculeC in moleculesC:

                        def generate_products_and_reactions(order):
                            """
                            order = (0, 1, 2) corresponds to reactants stored as A + B + C, etc.
                            """
                            _mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[order[0]])
                            _mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[order[1]])
                            _mappingsC = self.__matchReactantToTemplate(moleculeC, template.reactants[order[2]])

                            # Iterate over each pair of matches (A, B, C)
                            for _mapA in _mappingsA:
                                for _mapB in _mappingsB:
                                    for _mapC in _mappingsC:
                                        _reactantStructures = [moleculeA, moleculeB, moleculeC]
                                        _maps = [_mapA, _mapB, _mapC]
                                        # Reorder reactants in case we have a family with fewer reactant trees than
                                        # reactants and different reactant orders can produce different products
                                        _reactantStructures = [_reactantStructures[_i] for _i in order]
                                        _maps = [_maps[_i] for _i in order]
                                        try:
                                            _productStructures = self.__generateProductStructures(_reactantStructures,
                                                                                                  _maps,
                                                                                                  forward)
                                        except ForbiddenStructureException:
                                            pass
                                        else:
                                            if _productStructures is not None:
                                                _rxn = self.__createReaction(_reactantStructures,
                                                                             _productStructures,
                                                                             forward)
                                                if _rxn: rxnList.append(_rxn)

                        # Reactants stored as A + B + C
                        generate_products_and_reactions((0, 1, 2))

                        # Only check for swapped reactants if they are different
                        if reactants[1] is not reactants[2]:
                            # Reactants stored as A + C + B
                            generate_products_and_reactions((0, 2, 1))
                        if reactants[0] is not reactants[1]:
                            # Reactants stored as B + A + C
                            generate_products_and_reactions((1, 0, 2))
                        if reactants[0] is not reactants[2]:
                            # Reactants stored as C + B + A
                            generate_products_and_reactions((2, 1, 0))
                            if reactants[0] is not reactants[1] and reactants[1] is not reactants[2]:
                                # Reactants stored as C + A + B
                                generate_products_and_reactions((2, 0, 1))
                                # Reactants stored as B + C + A
                                generate_products_and_reactions((1, 2, 0))

        # If products is given, remove reactions from the reaction list that
        # don't generate the given products
        if products is not None:
            ensure_species(products, resonance=prod_resonance)

            rxnList0 = rxnList[:]
            rxnList = []
            for reaction in rxnList0:
            
                products0 = reaction.products[:] if forward else reaction.reactants[:]

                # For aromatics, generate aromatic resonance structures to accurately identify isomorphic species
                if prod_resonance:
                    for i, product in enumerate(products0):
                        if product.isCyclic:
                            aromaticStructs = generate_aromatic_resonance_structures(product)
                            if aromaticStructs:
                                products0[i] = aromaticStructs[0]

                # Skip reactions that don't match the given products
                if isomorphic_species_lists(products, products0):
                    rxnList.append(reaction)

        # Determine the reactant-product pairs to use for flux analysis
        # Also store the reaction template (useful so we can easily get the kinetics later)
        for reaction in rxnList:
            
            # Restore the labeled atoms long enough to generate some metadata
            for reactant in reaction.reactants:
                reactant.clearLabeledAtoms()
            for label, atom in reaction.labeledAtoms:
                atom.label = label
            
            # Generate metadata about the reaction that we will need later
            reaction.pairs = self.getReactionPairs(reaction)
            reaction.template = self.getReactionTemplateLabels(reaction)

            # Unlabel the atoms
            for label, atom in reaction.labeledAtoms:
                atom.label = ''
            
            # We're done with the labeled atoms, so delete the attribute
            del reaction.labeledAtoms

            # Mark reaction reversibility
            reaction.reversible = self.reversible
            
        # This reaction list has only checked for duplicates within itself, not
        # with the global list of reactions
        return rxnList

    def getReactionPairs(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, return the reactant-product pairs to use when
        performing flux analysis.
        """
        pairs = []
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
            elif reaction.reactants[1].containsLabeledAtom('*1'):
                if reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
        elif self.label.lower() in ['disproportionation', 'co_disproportionation', 'korcek_step1_cat']:
            # Hardcoding for disproportionation, co_disproportionation, korcek_step1_cat:
            # pair the reactant containing *1 with the product containing *1
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*1'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
            elif reaction.reactants[1].containsLabeledAtom('*1'):
                if reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
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
            elif reaction.reactants[1].containsLabeledAtom('*2'):
                if reaction.products[1].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
        elif self.label.lower() == 'baeyer-villiger_step1_cat':
            # Hardcoding for Baeyer-Villiger_step1_cat: pair the two reactants
            # with the Criegee intermediate and pair the catalyst with itself
            assert len(reaction.reactants) == 3 and len(reaction.products) == 2
            if reaction.reactants[0].containsLabeledAtom('*5'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                    pairs.append([reaction.reactants[2],reaction.products[0]])
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[0]])
            elif reaction.reactants[1].containsLabeledAtom('*5'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[2].containsLabeledAtom('*5'):
                if reaction.products[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                elif reaction.products[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
        elif self.label.lower() == 'baeyer-villiger_step2_cat':
            # Hardcoding for Baeyer-Villiger_step2_cat: pair the Criegee
            # intermediate with the two products and the catalyst with itself
            assert len(reaction.reactants) == 2 and len(reaction.products) == 3
            if reaction.products[0].containsLabeledAtom('*7'):
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                    pairs.append([reaction.reactants[0], reaction.products[0]])
            elif reaction.products[1].containsLabeledAtom('*7'):
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
            elif reaction.products[2].containsLabeledAtom('*7'):
                if reaction.reactants[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                elif reaction.reactants[1].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])

        if not pairs:
            logging.debug('Preset mapping missing for determining reaction pairs for family {0!s}, falling back to Reaction.generatePairs'.format(self.label))

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
        to use: 'group additivity' (new possible RMG-Py behavior) and 'rate rules' (old
        RMG-Java behavior, and default RMG-Py behavior).
        
        Returns a tuple (kinetics, entry):
        If it's estimated via 'rate rules' and an exact match is found in the tree,
        then the entry is returned as the second element of the tuple.
        But if an average is used, or the 'group additivity' method, then the tuple
        returned is (kinetics, None).
        """
        if method.lower() == 'group additivity':
            return self.estimateKineticsUsingGroupAdditivity(template, degeneracy), None
        elif method.lower() == 'rate rules':
            return self.estimateKineticsUsingRateRules(template, degeneracy)  # This returns kinetics and entry data
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
        entries = depository.entries.values()
        for entry in entries:
            if entry.item.isIsomorphic(reaction):
                kineticsList.append([deepcopy(entry.data), entry, entry.item.isIsomorphic(reaction, eitherDirection=False)])
        for kinetics, entry, is_forward in kineticsList:
            if kinetics is not None:
                kinetics.comment += "Matched reaction {0} {1} in {2}\nThis reaction matched rate rule {3}".format(entry.index, 
                                                      entry.label, 
                                                      depository.label,
                                                      '[{0}]'.format(';'.join([g.label for g in template])))
                kinetics.comment += "\nfamily: {}".format(self.label)
        return kineticsList
    
    def __selectBestKinetics(self, kineticsList):
        """
        For a given set of kinetics `kineticsList`, return the kinetics deemed
        to be the "best". This is determined to be the one with the lowest
        non-zero rank that occurs first (has the lowest index).
        """
        if any([x[1].rank == 0 for x in kineticsList]) and not all([x[1].rank == 0 for x in kineticsList]):
            kineticsList = [x for x in kineticsList if x[1].rank != 0]
        kineticsList.sort(key=lambda x: (x[1].rank, x[1].index))
        return kineticsList[0]
        
    def getKinetics(self, reaction, templateLabels, degeneracy=1, estimator='', returnAllKinetics=True):
        """
        Return the kinetics for the given `reaction` by searching the various
        depositories as well as generating a result using the user-specified `estimator`
        of either 'group additivity' or 'rate rules'.  Unlike
        the regular :meth:`getKinetics()` method, this returns a list of
        results, with each result comprising of

        1. the kinetics
        2. the source - this will be `None` if from a template estimate
        3. the entry  - this will be `None` if from a template estimate
        4. is_forward a boolean denoting whether the matched entry is in the same
           direction as the inputted reaction. This will always be True if using
           rates rules or group additivity. This can be `True` or `False` if using
           a depository

        If returnAllKinetics==False, only the first (best?) matching kinetics is returned.
        """
        kineticsList = []
        
        depositories = self.depositories[:]

        template = self.retrieveTemplate(templateLabels)
        
        # Check the various depositories for kinetics
        for depository in depositories:
            kineticsList0 = self.getKineticsFromDepository(depository, reaction, template, degeneracy)
            if len(kineticsList0) > 0 and not returnAllKinetics:
                kinetics, entry, is_forward = self.__selectBestKinetics(kineticsList0)
                return kinetics, depository, entry, is_forward
            else:
                for kinetics, entry, is_forward in kineticsList0:
                    kineticsList.append([kinetics, depository, entry, is_forward])
        
        # If estimator type of rate rules or group additivity is given, retrieve the kinetics. 
        if estimator:
            try:
                kinetics, entry = self.getKineticsForTemplate(template, degeneracy, method=estimator)
            except Exception:
                logging.error("Error getting kinetics for reaction {0!s}.\n{0!r}".format(reaction))
                raise

            if kinetics:
                if not returnAllKinetics:
                    return kinetics, estimator, entry, True
                kineticsList.append([kinetics, estimator, entry, True])
        # If no estimation method was given, prioritize rate rule estimation. 
        # If returning all kinetics, add estimations from both rate rules and group additivity.
        else:
            try:
                kinetics, entry = self.getKineticsForTemplate(template, degeneracy, method='rate rules')
                if not returnAllKinetics:
                    return kinetics, 'rate rules', entry, True
                kineticsList.append([kinetics, 'rate rules', entry, True])
            except KineticsError:
                # If kinetics were undeterminable for rate rules estimation, do nothing.
                pass
            
            try:
                kinetics2, entry2 = self.getKineticsForTemplate(template, degeneracy, method='group additivity')
                if not returnAllKinetics:
                    return kinetics, 'group additivity', entry2, True
                kineticsList.append([kinetics2, 'group additivity', entry2, True])
            except KineticsError:                
                # If kinetics were undeterminable for group additivity estimation, do nothing.
                pass
        
        if not returnAllKinetics:
            raise UndeterminableKineticsError(reaction)
        
        return kineticsList
    
    def estimateKineticsUsingGroupAdditivity(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using group additivity.
        
        Returns just the kinetics, or None.
        """
        warnings.warn("Group additivity is no longer supported and may be"
                      " removed in version 2.3.", DeprecationWarning)
        # Start with the generic kinetics of the top-level nodes
        kinetics = None
        root = self.getRootTemplate()
        kinetics = self.getKineticsForTemplate(root)
        
        if kinetics is None:
            #raise UndeterminableKineticsError('Cannot determine group additivity kinetics estimate for template "{0}".'.format(','.join([e.label for e in template])))
            return None
        else:
            kinetics = kinetics[0]
            
        # Now add in more specific corrections if possible
        return self.groups.estimateKineticsUsingGroupAdditivity(template, kinetics, degeneracy)        
        
    def estimateKineticsUsingRateRules(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.
        
        Returns a tuple (kinetics, entry) where `entry` is the database
        entry used to determine the kinetics only if it is an exact match,
        and is None if some averaging or use of a parent node took place.
        """    
        kinetics, entry  = self.rules.estimateKinetics(template, degeneracy)
                
        return kinetics, entry


    def getReactionTemplateLabels(self, reaction):
        """
        Retrieve the template for the reaction and 
        return the corresponding labels for each of the 
        groups in the template.
        """
        template = self.getReactionTemplate(reaction)
        
        templateLabels = []
        for entry in template:
            templateLabels.append(entry.label)

        return templateLabels

    def retrieveTemplate(self, templateLabels):
        """
        Reconstruct the groups associated with the 
        labels of the reaction template and 
        return a list.
        """
        template = []
        for label in templateLabels:
            template.append(self.groups.entries[label])

        return template

    def getLabeledReactantsAndProducts(self, reactants, products):
        """
        Given `reactants`, a list of :class:`Molecule` objects, and products, a list of 
        :class:`Molecule` objects, return two new lists of :class:`Molecule` objects with 
        atoms labeled: one for reactants, one for products. Returned molecules are totally 
        new entities in memory so input molecules `reactants` and `products` won't be affected.
        If RMG cannot find appropriate labels, (None, None) will be returned.
        """
        template = self.forwardTemplate
        reactants0 = [reactant.copy(deep=True) for reactant in reactants]

        if len(reactants0) == 1:
            molecule = reactants0[0]
            mappings = self.__matchReactantToTemplate(molecule, template.reactants[0])
            mappings = [[map0] for map0 in mappings]
            num_mappings = len(mappings)
            reactant_structures = [molecule]
        elif len(reactants0) == 2:
            moleculeA = reactants0[0]
            moleculeB = reactants0[1]
            # get mappings in forward direction
            mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[0])
            mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[1])
            mappings = list(itertools.product(mappingsA, mappingsB))
            # get mappings in the reverse direction
            mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[1])
            mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[0])
            mappings.extend(list(itertools.product(mappingsA, mappingsB)))

            reactant_structures = [moleculeA, moleculeB]
            num_mappings = len(mappingsA) * len(mappingsB)
        elif len(reactants0) == 3:
            moleculeA = reactants0[0]
            moleculeB = reactants0[1]
            moleculeC = reactants0[2]
            # Get mappings for all permutations of reactants
            mappings = []
            for order in itertools.permutations(range(3), 3):
                mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[order[0]])
                mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[order[1]])
                mappingsC = self.__matchReactantToTemplate(moleculeC, template.reactants[order[2]])
                mappings.extend(list(itertools.product(mappingsA, mappingsB, mappingsC)))

            reactant_structures = [moleculeA, moleculeB, moleculeC]
            num_mappings = len(mappingsA)*len(mappingsB)*len(mappingsC)
        else:
            raise IndexError('You have {0} reactants, which is unexpected!'.format(len(reactants)))

        for mapping in mappings:
            try:
                product_structures = self.__generateProductStructures(reactant_structures, mapping, forward=True)
            except ForbiddenStructureException:
                pass
            else:
                if product_structures is not None:
                    if isomorphic_species_lists(list(products), list(product_structures)):
                        return reactant_structures, product_structures
                    else:
                        continue

        # if there're some mapping available but cannot match the provided products
        # raise exception
        if num_mappings > 0:
            raise ActionError('Something wrong with products that RMG cannot find a match!')

        return None, None

    def addAtomLabelsForReaction(self, reaction, output_with_resonance = True):
        """
        Apply atom labels on a reaction using the appropriate atom labels from
        this reaction family.

        The reaction is modified in place containing species objects with the
        atoms labeled. If output_with_resonance is True, all resonance structures
        are generated with labels. If false, only the first resonance structure
        sucessfully able to map to the reaction is used. None is returned.
        """
        # make sure we start with reaction with species objects
        reaction.ensure_species(reactant_resonance=False, product_resonance=False)

        reactants = reaction.reactants
        products = reaction.products
        # ensure all species are independent references
        if len(reactants + products) > len(set([id(s) for s in reactants + products])):
            logging.debug('Copying reactants and products for reaction {} since they have identical species references'.format(reaction))
            # not all species are independent
            reactants = [s.copy(deep=True) for s in reactants]
            products = [s.copy(deep=True) for s in products]

        # get all possible pairs of resonance structures
        reactant_pairs = list(itertools.product(*[s.molecule for s in reaction.reactants]))
        product_pairs = list(itertools.product(*[s.molecule for s in reaction.products]))

        labeled_reactants, labeled_products = None, None
        # go through each combination of possible pairs
        for reactant_pair, product_pair in itertools.product(reactant_pairs, product_pairs):
            try:
                # see if we obtain proper labeling
                labeled_reactants, labeled_products = self.getLabeledReactantsAndProducts(reactant_pair, product_pair)
                if labeled_reactants is not None:
                    break
            except ActionError:
                # must have gotten the wrong pair
                pass
        if labeled_reactants is None or labeled_products is None:
            raise ActionError("Could not find labeled reactants for reaction {} from family {}.".format(reaction,self.label))

        # place the molecules in reaction's species object
        # this prevents overwriting of attributes of species objects by this method
        for index, species in enumerate(products):
            for labeled_molecule in labeled_products:
                if species.isIsomorphic(labeled_molecule):
                    species.molecule = [labeled_molecule]
                    reaction.products[index] = species
                    break
            else:
                raise ActionError('Could not find isomorphic molecule to fit the original product {} from reaction {}'.format(species, reaction))
        for index, species in enumerate(reactants):
            for labeled_molecule in labeled_reactants:
                if species.isIsomorphic(labeled_molecule):
                    species.molecule = [labeled_molecule]
                    reaction.reactants[index] = species
                    break
            else:
                raise ActionError('Could not find isomorphic molecule to fit the original reactant {} from reaction {}'.format(species, reaction))

        if output_with_resonance:
            # convert the molecules to species objects with resonance structures
            for species in reaction.reactants + reaction.products:
                species.generate_resonance_structures()

    def getTrainingDepository(self):
        """
        Returns the `training` depository from self.depositories
        """
        for depository in self.depositories:
            if depository.label.endswith('training'):
                return depository
        else:
            raise DatabaseError('Could not find training depository in family {0}.'.format(self.label))

            
    def addEntry(self,parent,grp,name):
        """
        Adds a group entry with parent parent
        group structure grp
        and group name name
        """
        ind = len(self.groups.entries)-1
        entry = Entry(index=ind,label=name,item=grp,parent=parent)
        self.groups.entries[name] = entry
        self.rules.entries[name] = []
        if entry.parent:
            entry.parent.children.append(entry)
    
    def getEntriesReactions(self,template):
        """
        retrieves all training reactions whose kinetics
        are associated with the entry template
        """
        entries = self.rules.entries[template]
        tentries = self.getTrainingDepository().entries
        rxnEntries = []
        for entry in entries:
            rxnEntries.append(tentries[entry.index].item)
        return rxnEntries
    
    def getTemplateKinetics(self,template):
        """
        retrives a list of all the kinetics objects
        associated with a given template
        """
        return [entry.data for entry in self.rules.entries[template]]
    
    def splitReactions(self,rxns,oldlabel,newgrp):
        """
        divides the reactions in rxns between the new
        group structure newgrp and the old structure with 
        label oldlabel
        returns a list of reactions associated with the new group
        the list of reactions associated with the old group
        and a list of the indices of all of the reactions
        associated with the new group
        """
        new = []
        comp = []
        newInds = []
        kinetics = self.getTemplateKinetics(oldlabel)
        
        for i,rxn in enumerate(rxns):
            reactants = rxn.reactants if rxn.is_forward else rxn.products
            rmol = reactants[0].molecule[0]
            for r in reactants[1:]:
                rmol.merge(r.molecule[0])
            
            if rmol.isSubgraphIsomorphic(newgrp,generateInitialMap=True, saveOrder=True):
                new.append(kinetics[i])
                newInds.append(i)
            else:
                comp.append(kinetics[i])
        
        return new,comp,newInds    
    
    def evalExt(self,parent,ext,extname,obj=None,T=1000.0):
        """
        evaluates the objective function obj
        for the extension ext with name extname to the parent entry parent
        """
        rxns = self.getEntriesReactions(parent.label)
        new,old,newInds = self.splitReactions(rxns,parent.label,ext)
        if len(new) == 0:
            return np.inf,False
        elif len(old) == 0:
            return np.inf,True
        else:
            if obj:
                ob,boo = getObjectiveFunction(new,old,obj,T=T)
            else:
                ob,boo = getObjectiveFunction(new,old,T=T)
            return ob,True
    
    def getExtensionEdge(self,parent,obj,T):
        """
        finds the set of all extension groups to parent such that
        1) the extension group divides the set of reactions under parent
        2) No generalization of the extension group divides the set of reactions under parent
                    
        We find this by generating all possible extensions of the initial group.  Extensions that split reactions are added
        to the list.  All extensions that do not split reactions and do not create bonds are ignored 
        (although those that match every reaction are labeled so we don't search them twice).  Those that match
        all reactions and involve bond creation undergo this process again.  
        
        Principle:  Say you have two elementary changes to a group ext1 and ext2 if applying ext1 and ext2 results in a 
        split at least one of ext1 and ext2 must result in a split
        
        Speed of this algorithm relies heavily on searching non bond creation dimensions once.
        """
        outExts = [[]]
        grps = [parent.item]
        names = [parent.label]
        atmInds = [None]
        
        while grps != []:
            grp = grps[-1]
            
            if atmInds[-1]:
                if len(atmInds[-1]) == 1:
                    exts = grp.getExtensions(basename=names[-1],atmInd=atmInds[-1][0])
                elif len(atmInds[-1]) == 2:
                    exts = grp.getExtensions(basename=names[-1],atmInd=atmInds[-1][0],atmInd2=atmInds[-1][1])
            else:
                exts = grp.getExtensions(basename=names[-1])
            
            regDict = dict()
            extInds = []
            for i,(grp2,grpc,name,typ,indc) in enumerate(exts):

                if typ != 'intNewBondExt' and typ != 'extNewBondExt' and (typ,indc) not in regDict.keys():
                    regDict[(typ,indc)] = ([],[])
                val,boo = self.evalExt(parent,grp2,name,obj,T)
                    
                if val != np.inf:
                    outExts[-1].append(exts[i]) #this extension splits reactions (optimization dim)
                    if typ == 'atomExt':
                        regDict[(typ,indc)][0].extend(grp2.atoms[indc[0]].atomType)
                    elif typ == 'elExt':
                        regDict[(typ,indc)][0].extend(grp2.atoms[indc[0]].radicalElectrons)
                    elif typ == 'bondExt':
                        regDict[(typ,indc)][0].extend(grp2.getBond(grp2.atoms[indc[0]],grp2.atoms[indc[1]]).order)
                        
                elif boo: #this extension matches all reactions (regularization dim)
                    if typ == 'intNewBondExt' or typ == 'extNewBondExt':
                        extInds.append(i)  #these are bond formation extensions, we want to expand these until we get splits 
                    elif typ == 'atomExt':
                        regDict[(typ,indc)][0].extend(grp2.atoms[indc[0]].atomType)
                    elif typ == 'elExt':
                        regDict[(typ,indc)][0].extend(grp2.atoms[indc[0]].radicalElectrons)
                    elif typ == 'bondExt':
                        regDict[(typ,indc)][0].extend(grp2.getBond(grp2.atoms[indc[0]],grp2.atoms[indc[1]]).order)

                else:                    
                    #this extension matches no reactions
                    if typ == 'atomExt':
                        regDict[(typ,indc)][1].extend(grp2.atoms[indc[0]].atomType)
                    elif typ == 'elExt':
                        regDict[(typ,indc)][1].extend(grp2.atoms[indc[0]].radicalElectrons)
                    elif typ == 'bondExt':
                        regDict[(typ,indc)][1].extend(grp2.getBond(grp2.atoms[indc[0]],grp2.atoms[indc[1]]).order)
            
                    
            for typr,indcr in regDict.keys(): #have to label the regularization dimensions in all relevant groups
                regVal = regDict[(typr,indcr)][0]
                #parent
                if typr != 'intNewBondExt' and typr != 'extNewBondExt': #these dimensions should be regularized
                    if typr == 'atomExt':
                        grp.atoms[indcr[0]].reg_dim_atm = regVal
                    elif typr == 'elExt':
                        grp.atoms[indcr[0]].reg_dim_u = regVal
                    elif typr == 'bondExt':
                        atms = grp.atoms
                        bd = grp.getBond(atms[indcr[0]],atms[indcr[1]])
                        bd.reg_dim = regVal
                            
                #extensions being sent out
                if typr != 'intNewBondExt' and typr != 'extNewBondExt': #these dimensions should be regularized
                    for grp2,grpc,name,typ,indc in outExts[-1]: #returned groups
                        if typr == 'atomExt':
                            grp2.atoms[indcr[0]].reg_dim_atm = regVal
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_atm = regVal
                        elif typr == 'elExt':
                            grp2.atoms[indcr[0]].reg_dim_u = regVal
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_u = regVal
                        elif typr == 'bondExt':
                            atms = grp2.atoms
                            bd = grp2.getBond(atms[indcr[0]],atms[indcr[1]])
                            bd.reg_dim = regVal
                            if grpc:
                                atms = grpc.atoms
                                bd = grp2.getBond(atms[indcr[0]],atms[indcr[1]])
                                bd.reg_dim = regVal
            
            #extensions being expanded
            for typr,indcr in regDict.keys(): #have to label the regularization dimensions in all relevant groups
                regVal = regDict[(typr,indcr)][0]
                if typr != 'intNewBondExt' and typr != 'extNewBondExt': #these dimensions should be regularized
                    for ind2 in extInds: #groups for expansion
                        grp2,grpc,name,typ,indc = exts[ind2]
                        if typr == 'atomExt':
                            grp2.atoms[indcr[0]].reg_dim_atm = regVal
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_atm = regVal
                        elif typr == 'elExt':
                            grp2.atoms[indcr[0]].reg_dim_u = regVal
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_u = regVal
                        elif typr == 'bondExt':
                            atms = grp2.atoms
                            bd = grp2.getBond(atms[indcr[0]],atms[indcr[1]])
                            bd.reg_dim = regVal
                            if grpc:
                                atms = grpc.atoms
                                bd = grp2.getBond(atms[indcr[0]],atms[indcr[1]])
                                bd.reg_dim = regVal
            
            outExts.append([])
            grps.pop()
            names.pop()
            atmInds.pop()
            
            for ind in extInds: #collect the groups to be expanded
                grpr,grpcr,namer,typr,indcr = exts[ind]
                grps.append(grpr)
                names.append(namer)
                atmInds.append(indcr)
        
        out = []
        for x in outExts: #compile all of the valid extensions together, may be some duplicates here, but I don't think it's currently worth identifying them
            out.extend(x)
        
        return out
            
                    
    def extendNode(self,parent,thermoDatabase=None,obj=None,T=1000.0):
        """
        Constructs an extension to the group parent based on evaluation 
        of the objective function obj
        """
        
        exts = self.getExtensionEdge(parent,obj=obj,T=T)
        
        vals = []
        for grp,grpc,name,typ,einds in exts:
            val,boo = self.evalExt(parent,grp,name,obj,T)
            vals.append(val) 
            
        min_val = min(vals)
        
        min_ind = vals.index(min_val)
        
        ext = exts[min_ind]
        
        extname = ext[2]
        
        self.addEntry(parent,ext[0],extname)
        
        complement = not ext[1] is None
        
        if complement:
            frags = extname.split('_')
            frags[-1] = 'N-'+frags[-1]
            cextname = ''
            for k in frags:
                cextname += k
                cextname += '_'
            cextname = cextname[:-1]
    
            self.addEntry(parent,ext[1],cextname)
        
        rxns = self.getEntriesReactions(parent.label)
        new,left,newInds = self.splitReactions(rxns,parent.label,ext[0])
        
        compEntries = []
        newEntries = []

        for i,entry in enumerate(self.rules.entries[parent.label]):
            if i in newInds:
                entry.label = extname
                newEntries.append(entry)
            else:
                if complement:
                    entry.label = cextname
                compEntries.append(entry)
        
        self.rules.entries[extname] = newEntries
        
        if complement:
            self.rules.entries[parent.label] = []
            self.rules.entries[cextname] = compEntries
        else:
            self.rules.entries[parent.label] = compEntries
            
        return
    
    def generateTree(self,obj=None,thermoDatabase=None,T=1000.0):
        """
        Generate a tree by greedy optimization based on the objective function obj
        the optimization is done by iterating through every group and if the group has
        more than one training reaction associated with it a set of potential more specific extensions 
        are generated and the extension that optimizing the objective function combination is chosen 
        and the iteration starts over at the beginning
        
        additionally the tree structure is simplified on the fly by removing groups that have no kinetics data associated
        if their parent has no kinetics data associated and they either have only one child or
        have two children one of which has no kinetics data and no children
        (its parent becomes the parent of its only relevant child node)
        """
        boo = True #if the for loop doesn't break becomes false and the while loop terminates
        while boo:
            for entry in self.groups.entries.itervalues():
                if not isinstance(entry.item, Group): #skip logic nodes
                    continue
                if entry.index != -1 and len(self.rules.entries[entry.label])>1:
                    self.extendNode(entry,thermoDatabase,obj,T)
                    break
                elif entry.parent is None or entry.parent.parent is None or not isinstance(entry.parent.item,Group) or not isinstance(entry.parent.parent.item,Group):
                    pass
                elif len(self.rules.entries[entry.parent.label])>0 or len(self.rules.entries[entry.label])>0:
                    pass
                elif len(entry.parent.children) == 1:
                    label = entry.parent.label
                    entry.parent.parent.children.remove(entry.parent)
                    entry.parent = entry.parent.parent
                    entry.parent.children.append(entry)
                    del self.groups.entries[label]   
                    del self.rules.entries[label]
                    break
                elif len(entry.parent.children) == 2: 
                    child = [c for c in entry.parent.children if c != entry][0]
                    if len(self.rules.entries[child.label]) == 0 and len(child.children) == 0:
                        label = entry.parent.label
                        entry.parent.parent.children.remove(entry.parent)
                        entry.parent.parent.children.append(entry)
                        entry.parent = entry.parent.parent
                        del self.groups.entries[label]   
                        del self.rules.entries[label]
                        clabel = child.label
                        del self.groups.entries[clabel]
                        del self.rules.entries[clabel]
                        break
            else:
                boo = False
            
            #fix indicies
            iters = 0
            for entry in self.groups.entries.itervalues():
                if entry.index != -1:
                    entry.index = iters
                    iters += 1
        
        return
    
    def simpleRegularization(self, node):
        """
        Simplest regularization algorithm
        All nodes are made as specific as their descendant reactions
        Training reactions are assumed to not generalize 
        For example if an particular atom at a node is Oxygen for all of its
        descendent reactions a reaction where it is Sulfur will never hit that node
        unless it is the top node even if the tree did not split on the identity 
        of that atom
        """
        grp = node.item
        
        if isinstance(node.item,Group):
            for i,atm1 in enumerate(grp.atoms):
                if atm1.reg_dim_atm != [] and set(atm1.reg_dim_atm) != set(atm1.atomType):
                    self.extendRegularization(node,[i],atm1.reg_dim_atm,'atomtype')
                if atm1.reg_dim_u != [] and set(atm1.reg_dim_u) != set(atm1.radicalElectrons):
                    self.extendRegularization(node,[i],atm1.reg_dim_u,'unpaired')
                for j,atm2 in enumerate(grp.atoms[i:]):
                    if grp.hasBond(atm1,atm2):
                        bd = grp.getBond(atm1,atm2)
                        if bd.reg_dim != [] and set(bd.reg_dim) != set(bd.order):
                            self.extendRegularization(node,[i,j],bd.reg_dim,'bond')    
        
        for child in node.children:
            self.simpleRegularization(child)
            
    def extendRegularization(self, node, inds, regs, typ):
        """
        Applies a regularization down the tree from a given parent node
        """
        grp = node.item
        
        if isinstance(grp,Group):
            if typ == 'atomtype':
                grp.atoms[inds[0]].atomType = list(set(grp.atoms[inds[0]].atomType) & set(regs))
                for child in node.children:
                    self.extendRegularization(child,inds,regs,typ)
            elif typ == 'unpaired':
                grp.atoms[inds[0]].radicalElectrons = list(set(grp.atoms[inds[0]].radicalElectrons) & set(regs))
                for child in node.children:
                    self.extendRegularization(child,inds,regs,typ)
            elif typ == 'bond':
                bd = grp.getBond(grp.atoms[inds[0]],grp.atoms[inds[1]])
                bd.order = list(set(bd.order) & set(regs))
                for child in node.children:
                    self.extendRegularization(child,inds,regs,typ)
            else:
                raise ValueError('regularization type of {0} is unimplemented'.format(typ))
    
    def regularize(self, regularization=simpleRegularization):
        """
        Regularizes the tree according to the regularization function regularization
        """
        regularization(self,self.getRootTemplate()[0])
        
    def prepareTreeForGeneration(self,thermoDatabase=None):
        """
        clears groups and rules in the tree, generates an appropriate
        root group to start from and then reads training reactions
        Note this only works if a single top node (not a logic node)
        can be generated
        """
        #find the starting node
        grp = None
        
        rtmps = self.getRootTemplate()
        
        if not isinstance(rtmps[0].item,Group):
            raise ValueError('each tree top node must be a group not a logic node to prepare the tree automatically')
        
        for ent in rtmps:
            if grp is None:
                grp = ent.item
            else:
                grp.mergeGroups(ent.item)
        
        
        #clear everything
        self.groups.entries = {x.label:x for x in self.groups.entries.itervalues() if x.index == -1}
        self.rules.entries = OrderedDict()
        
        #add the starting node
        self.addEntry(None,grp,'Root')
        self.groups.top = [self.groups.entries['Root']]
        self.forwardTemplate.reactants = [self.groups.entries['Root']]

        #fill with training reactions
        self.addKineticsRulesFromTrainingSet(thermoDatabase)
        
        return
    
    def saveGeneratedTree(self,path=None):
        """
        clears the rules and saves the family to its 
        current location in database
        """
        if path is None:
            path = settings['database.directory']
            path = os.path.join(path,'kinetics','families',self.label)
        
        self.rules.entries = OrderedDict() #have to clear the new rules made
        
        self.save(path)
    
    def retrieveOriginalEntry(self, templateLabel):
        """
        Retrieves the original entry, be it a rule or training reaction, given
        the template label in the form 'group1;group2' or 'group1;group2;group3'
        
        Returns tuple in the form
        (RateRuleEntry, TrainingReactionEntry)
        
        Where the TrainingReactionEntry is only present if it comes from a training reaction
        """
        templateLabels = templateLabel.split()[-1].split(';')
        template = self.retrieveTemplate(templateLabels)
        rule = self.getRateRule(template)
        if 'From training reaction' in rule.data.comment:
            trainingIndex = int(rule.data.comment.split()[3])
            trainingDepository = self.getTrainingDepository()
            return rule, trainingDepository.entries[trainingIndex]
        else:
            return rule, None
        
    def getSourcesForTemplate(self, template):
        """
        Returns the set of rate rules and training reactions used to average this `template`.  Note that the tree must be
        averaged with verbose=True for this to work.
        
        Returns a tuple of
        rules, training
        
        where rules are a list of tuples containing 
        the [(original_entry, weight_used_in_average), ... ]
        
        and training is a list of tuples containing
        the [(rate_rule_entry, training_reaction_entry, weight_used_in_average),...]
        """
        import re

        def assignWeightsToEntries(entryNestedList, weightedEntries, N = 1):
            """
            Assign weights to an average of average nested list. Where N is the 
            number of values being averaged recursively.  
            """
            N = len(entryNestedList)*N
            for entry in entryNestedList:
                if isinstance(entry, list):
                    assignWeightsToEntries(entry, weightedEntries, N)
                else:
                    weightedEntries.append((entry,1/float(N)))
            return weightedEntries
        
        
        kinetics, entry = self.estimateKineticsUsingRateRules(template)
        if entry:
            return [(entry,1)], []   # Must be a rate rule 
        else:
            # The template was estimated using an average or another node
            rules = []
            training = []
            
            lines = kinetics.comment.split('\n')
            
            lines = [line for line in lines if not line.startswith('Euclid') and not line.startswith('family:')] #remove the Euclidean distance and family lines to help parser
            
            # Discard the last line, unless it's the only line!
            # The last line is 'Estimated using ... for rate rule (originalTemplate)'
            #if from training reaction is in the first line append it to the end of the second line and skip the first line
            if not 'Average of' in kinetics.comment:
                if 'From training reaction' in lines[0]:
                    comment = lines[1]
                else:
                    comment = lines[0]
                if comment.startswith('Estimated using template'):
                    tokenTemplateLabel = comment.split()[3][1:-1]
                    ruleEntry, trainingEntry = self.retrieveOriginalEntry(tokenTemplateLabel) 
                    if trainingEntry:
                        training.append((ruleEntry,trainingEntry,1))   # Weight is 1
                    else:
                        rules.append((ruleEntry,1))
                else:
                    raise ValueError('Could not parse unexpected line found in kinetics comment: {}'.format(comment))
            else:
                comment = ' '.join(lines[:-1])
                # Clean up line for exec
                evalCommentString = re.sub(r" \+ ", ",",                        # any remaining + signs
                                    re.sub(r"Average of ", "",                  # average of averages
                                    re.sub(r"Average of \[(?!Average)", "['",   # average of groups
                                    re.sub(r"(\w|\))]", r"\1']",                # initial closing bracket
                                    re.sub(r"(?<=[\w)]) \+ (?=Average)", "',",  # + sign between non-average and average
                                    re.sub(r"(?<=]) \+ (?!Average)", ",'",      # + sign between average and non-average
                                    re.sub(r"(?<!]) \+ (?!Average)", "','",     # + sign between non-averages
                                    comment)))))))

                entryNestedList = eval(evalCommentString)
                
                weightedEntries = assignWeightsToEntries(entryNestedList, [])
                
                
                rules = {}
                training = {}
                
                for tokenTemplateLabel, weight in weightedEntries:
                    if 'From training reaction' in tokenTemplateLabel:
                        tokenTemplateLabel = tokenTemplateLabel.split()[-1]
                    ruleEntry, trainingEntry = self.retrieveOriginalEntry(tokenTemplateLabel)
                    if trainingEntry:
                        if (ruleEntry, trainingEntry) in training:
                            training[(ruleEntry, trainingEntry)] += weight
                        else:
                            training[(ruleEntry, trainingEntry)] = weight
                    else:
                        if ruleEntry in rules:
                            rules[ruleEntry] += weight
                        else:
                            rules[ruleEntry] = weight
                # Each entry should now only appear once    
                training = [(k[0],k[1],v) for k,v in training.items()]
                rules = rules.items()
                
            return rules, training

    def extractSourceFromComments(self, reaction):
        """
        Returns the rate rule associated with the kinetics of a reaction by parsing the comments.
        Will return the template associated with the matched rate rule.
        Returns a tuple containing (Boolean_Is_Kinetics_From_Training_reaction, Source_Data)
        
        For a training reaction, the Source_Data returns::

            [Family_Label, Training_Reaction_Entry, Kinetics_In_Reverse?]
        
        For a reaction from rate rules, the Source_Data is a tuple containing::

            [Family_Label, {'template': originalTemplate,
                            'degeneracy': degeneracy,
                            'exact': boolean_exact?,
                            'rules': a list of (original rate rule entry, weight in average)
                            'training': a list of (original rate rule entry associated with training entry, original training entry, weight in average)}]


        where Exact is a boolean of whether the rate is an exact match, Template is
        the reaction template used, RateRules is a list of the rate rule entries containing
        the kinetics used, and TrainingReactions are ones that have created rules used in the estimate.
        """
        import re
        lines = reaction.kinetics.comment.split('\n')

        exact = False
        template = None
        rules = None
        trainingEntries = None
        degeneracy = 1

        regex = "\[(.*)\]" # only hit outermost brackets
        for line in lines:
            if line.startswith('Matched'):
                # Source of the kinetics is from training reaction
                trainingReactionIndex = int(line.split()[2])
                depository  = self.getTrainingDepository()
                trainingEntry = depository.entries[trainingReactionIndex]
                # Perform sanity check that the training reaction's label matches that of the comments
                if trainingEntry.label not in line:
                    raise AssertionError('Reaction {0} uses kinetics from training reaction {1} but does not match the training reaction {1} from the {2} family.'.format(reaction,trainingReactionIndex,self.label))
                
                # Sometimes the matched kinetics could be in the reverse direction..... 
                if reaction.isIsomorphic(trainingEntry.item, eitherDirection=False):
                    reverse=False
                else:
                    reverse=True
                return True, [self.label, trainingEntry, reverse]

            elif line.startswith('Exact match'):
                exact = True
            elif line.startswith('Estimated'):
                pass
            elif line.startswith('Multiplied by'):
                degeneracy = float(line.split()[-1])

        # Extract the rate rule information 
        fullCommentString = reaction.kinetics.comment.replace('\n', ' ')
        
        # The rate rule string is right after the phrase 'for rate rule'
        rateRuleString = fullCommentString.split("for rate rule",1)[1].split()[0]
        
        if rateRuleString[0] == '[':
            templateLabel = re.split(regex, rateRuleString)[1]
        else:
            templateLabel = rateRuleString #if has the line 'From training reaction # for rate rule node1;node2'
            
        template = self.retrieveTemplate(templateLabel.split(';'))
        rules, trainingEntries = self.getSourcesForTemplate(template)
        

        if not template:
            raise ValueError('Could not extract kinetics source from comments for reaction {}.'.format(reaction))
        
        sourceDict = {'template':template, 'degeneracy':degeneracy, 'exact':exact, 
                       'rules':rules,'training':trainingEntries }

        # Source of the kinetics is from rate rules
        return False, [self.label, sourceDict]

    def getBackboneRoots(self):
        """
        Returns: the top level backbone node in a unimolecular family.
        """

        backboneRoots = [entry for entry in self.groups.top if entry in self.forwardTemplate.reactants]
        return backboneRoots

    def getEndRoots(self):
        """
        Returns: A list of top level end nodes in a unimolecular family
        """

        endRoots = [entry for entry in self.groups.top if entry not in self.forwardTemplate.reactants]
        return endRoots

    def getTopLevelGroups(self, root):
        """
        Returns a list of group nodes that are the highest in the tree starting at node "root".
        If "root" is a group node, then it will return a single-element list with "root".
        Otherwise, for every child of root, we descend until we find no nodes with logic
        nodes. We then return a list of all group nodes found along the way.
        """

        groupList = [root]
        allGroups = False

        while not allGroups:
            newGroupList = []
            for entry in groupList:
                if isinstance(entry.item,Group):
                    newGroupList.append(entry)
                else:
                    newGroupList.extend(entry.children)
            groupList = newGroupList
            allGroups = all([isinstance(entry.item, Group) for entry in groupList])

        return groupList

def informationGain(ks1,ks2):
    """
    calculates the information gain as the sum of the products of the standard deviations at each
    node and the number of reactions at that node
    """
    return len(ks1)*np.std(ks1)+len(ks2)*np.std(ks2)
 
def getObjectiveFunction(kinetics1,kinetics2,obj=informationGain,T=1000.0):
    """
    Returns the value of four potential objective functions to minimize
    Uncertainty = N1*std(Ln(k))_1 + N1*std(Ln(k))_1
    Mean difference: -abs(mean(Ln(k))_1-mean(Ln(k))_2)
    Error using mean: Err_1 + Err_2
    Split: abs(N1-N2)
    """
    ks1 = np.array([np.log(k.getRateCoefficient(T)) for k in kinetics1])
    ks2 = np.array([np.log(k.getRateCoefficient(T)) for k in kinetics2])
    N1 = len(ks1)
    
    return obj(ks1,ks2), N1 == 0

