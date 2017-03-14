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

import os.path
import logging
import codecs
from copy import deepcopy

from rmgpy.constraints import failsSpeciesConstraints
from rmgpy.data.base import Database, Entry, LogicNode, LogicOr, ForbiddenStructures,\
                            ForbiddenStructureException, getAllCombinations
from rmgpy.reaction import Reaction
from rmgpy.kinetics import Arrhenius, ArrheniusEP
from rmgpy.molecule import Bond, GroupBond, Group, Molecule, ActionError
from rmgpy.species import Species
from rmgpy.molecule.resonance import generateAromaticResonanceStructures

from .common import KineticsError, UndeterminableKineticsError, saveEntry
from .depository import KineticsDepository
from .groups import KineticsGroups
from .rules import KineticsRules

################################################################################

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
                reverse=None
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
        self.reverse = reverse

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
                                   self.estimator,
                                   self.reverse,
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
    `forwardTemplate`   :class:`Reaction`               The forward reaction template
    `forwardRecipe`     :class:`ReactionRecipe`         The steps to take when applying the forward reaction to a set of reactants
    `reverseTemplate`   :class:`Reaction`               The reverse reaction template
    `reverseRecipe`     :class:`ReactionRecipe`         The steps to take when applying the reverse reaction to a set of reactants
    `forbidden`         :class:`ForbiddenStructures`    (Optional) Forbidden product structures in either direction
    `ownReverse`        `Boolean`                       It's its own reverse?
    'boundaryAtoms'     list                            Labels which define the boundaries of end groups in backbone/end families
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
                 shortDesc='',
                 longDesc='',
                 forwardTemplate=None,
                 forwardRecipe=None,
                 reverseTemplate=None,
                 reverseRecipe=None,
                 forbidden=None,
                 boundaryAtoms = None
                 ):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
        self.reverse = reverse
        self.forwardTemplate = forwardTemplate
        self.forwardRecipe = forwardRecipe
        self.reverseTemplate = reverseTemplate
        self.reverseRecipe = reverseRecipe
        self.forbidden = forbidden
        self.ownReverse = forwardTemplate is not None and reverseTemplate is None
        self.boundaryAtoms = boundaryAtoms
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
            
        self.rules.saveOld(path, self)
            
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
        local_context['reverse'] = None
        local_context['boundaryAtoms'] = None

        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        logging.debug("Loading kinetics family groups from {0}".format(os.path.join(path, 'groups.py')))
        Database.load(self.groups, os.path.join(path, 'groups.py'), local_context, global_context)
        self.name = self.label
        self.boundaryAtoms = local_context.get('boundaryAtoms', None)

        # Generate the reverse template if necessary
        self.forwardTemplate.reactants = [self.groups.entries[label] for label in self.forwardTemplate.reactants]
        if self.ownReverse:
            self.forwardTemplate.products = self.forwardTemplate.reactants[:]
            self.reverseTemplate = None
            self.reverseRecipe = None
        else:
            self.reverse = local_context.get('reverse', None)
            if self.reverse is None:
                self.reverse = '{0}_reverse'.format(self.label)
            self.forwardTemplate.products = self.generateProductTemplate(self.forwardTemplate.reactants)
            self.reverseTemplate = Reaction(reactants=self.forwardTemplate.products, products=self.forwardTemplate.reactants)
            self.reverseRecipe = self.forwardRecipe.getReverse()
        
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
    
    def saveTrainingReactions(self, reactions, reference=None, referenceType='', shortDesc='', rank=3):
        """
        This function takes a list of reactions appends it to the training reactions file.  It ignores the existence of
        duplicate reactions.  
        
        The rank for each new reaction's kinetics is set to a default value of 3 unless the user specifies differently 
        for those reactions.
        
        For each entry, the long description is imported from the kinetics comment. 
        """ 
        from rmgpy import settings

        training_path = os.path.join(settings['database.directory'], 'kinetics', 'families', \
            self.label, 'training')

        directory_file = os.path.join(training_path, 'dictionary.txt')

        # Load the old set of the species of the training reactions
        speciesDict = Database().getSpecies(directory_file)

        # add new unique species with labeledAtoms into speciesDict
        for rxn in reactions:
            for spec in (rxn.reactants + rxn.products):
                for ex_spec_label in speciesDict:
                    ex_spec = speciesDict[ex_spec_label]
                    if ex_spec.molecule[0].getFormula() != spec.molecule[0].getFormula():
                        continue
                    else:
                        spec_labeledAtoms = spec.molecule[0].getLabeledAtoms()
                        ex_spec_labeledAtoms = ex_spec.molecule[0].getLabeledAtoms()
                        initialMap = {}
                        try:
                            for atomLabel in spec_labeledAtoms:
                                initialMap[spec_labeledAtoms[atomLabel]] = ex_spec_labeledAtoms[atomLabel]
                        except KeyError:
                            # atom labels did not match, therefore not a match
                            continue
                        if spec.molecule[0].isIsomorphic(ex_spec.molecule[0],initialMap):
                            spec.label = ex_spec.label
                            break
                else:# no isomorphic existing species found
                    spec_formula = spec.molecule[0].getFormula()
                    if spec_formula not in speciesDict:
                        spec.label = spec_formula
                    else:
                        index = 2
                        while (spec_formula + '-{}'.format(index)) in speciesDict:
                            index += 1
                        spec.label = spec_formula + '-{}'.format(index)
                    speciesDict[spec.label] = spec

        training_file = open(os.path.join(settings['database.directory'], 'kinetics', 'families', \
            self.label, 'training', 'reactions.py'), 'a')
        training_file.write("\n\n")

        # get max reaction entry index from the existing training data
        try:
            depository = self.getTrainingDepository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Starting a new one')
            depository = KineticsDepository()
            self.depositories.append(depository)
        
        trainingDatabase = depository
        indices = [entry.index for entry in trainingDatabase.entries.values()]
        if indices:
            maxIndex = max(indices)
        else:
            maxIndex = 0
        # save new reactions to reactions.py
        for i, reaction in enumerate(reactions):    
            longDesc = 'Taken from entry: {0}'.format(reaction.kinetics.comment)
            reaction.kinetics.comment = ''
            entry = Entry(
                index = maxIndex+i+1,
                label = str(reaction),
                item = reaction,
                data = reaction.kinetics,
                reference = reference,
                referenceType = referenceType,
                shortDesc = unicode(shortDesc),
                longDesc = unicode(longDesc),
                rank = rank,
                )
            
            self.saveEntry(training_file, entry)
        training_file.close()

        # save species to dictionary
        with open(directory_file, 'w') as f:
            for label in speciesDict.keys():
                f.write(speciesDict[label].molecule[0].toAdjacencyList(label=label, removeH=False))
                f.write('\n')
        f.close()

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
                        print struct.toAdjacencyList()
                        print s.toAdjacencyList()
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
            
            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            
            new_entry = Entry(
                index = index,
                label = ';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template],
                                                   products=[]),
                data = ArrheniusEP(
                    A = deepcopy(data.A),
                    n = deepcopy(data.n),
                    alpha = 0,
                    E0 = deepcopy(data.Ea),
                    Tmin = deepcopy(data.Tmin),
                    Tmax = deepcopy(data.Tmax),
                    comment = "{0} from training reaction {1}".format(';'.join([g.label for g in template]), entry.index),
                ),
                rank = entry.rank,
                reference=entry.reference,
                shortDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.shortDesc,
                longDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.longDesc,
            )
            new_entry.data.A.value_si /= entry.item.degeneracy
            try:
                self.rules.entries[new_entry.label].append(new_entry)
            except KeyError:
                self.rules.entries[new_entry.label] = [new_entry]
            index += 1
        
        # Process the entries that are stored in the reverse direction of the
        # family definition
        for entry in reverse_entries:
            
            assert isinstance(entry.data, Arrhenius)
            data = deepcopy(entry.data)
            data.changeT0(1)
            # Estimate the thermo for the reactants and products
            # trainingSet=True used later to does not allow species to match a liquid phase library and get corrected thermo which will affect reverse rate calculation
            item = Reaction(reactants=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label) for m in entry.item.reactants],
                             products=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label) for m in entry.item.products])
            for reactant in item.reactants:
                reactant.generateResonanceIsomers()
                reactant.thermo = thermoDatabase.getThermoData(reactant, trainingSet=True) 
            for product in item.products:
                product.generateResonanceIsomers()
                product.thermo = thermoDatabase.getThermoData(product,trainingSet=True)
            # Now that we have the thermo, we can get the reverse k(T)
            item.kinetics = data
            data = item.generateReverseRateCoefficient()
            
            item = Reaction(reactants=[m.molecule[0].copy(deep=True) for m in entry.item.products], products=[m.molecule[0].copy(deep=True) for m in entry.item.reactants])
            template = self.getReactionTemplate(item)
            item.degeneracy = self.calculateDegeneracy(item)
            
            new_entry = Entry(
                index = index,
                label = ';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template],
                                                   products=[]),
                data = ArrheniusEP(
                    A = deepcopy(data.A),
                    n = deepcopy(data.n),
                    alpha = 0,
                    E0 = deepcopy(data.Ea),
                    Tmin = deepcopy(data.Tmin),
                    Tmax = deepcopy(data.Tmax),
                    comment = "{0} from training reaction {1}".format(';'.join([g.label for g in template]), entry.index),
                ),
                rank = entry.rank,
                reference=entry.reference,
                shortDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.shortDesc,
                longDesc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.longDesc,
            )
            new_entry.data.A.value_si /= item.degeneracy
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

        # For molecules, we need to assign atom indices in order to track atoms for isomorphism checking later
        if isinstance(reactantStructure, Molecule):
            reactantStructure.assignAtomIDs()

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
                raise KineticsError('Unable to change labels from "*" to "*1" and "*2" for reaction family {0}.'.format(label))

        # Generate the product structure by applying the recipe
        if forward:
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
                        
            elif label == 'H_shift_cyclopentadiene':
                # Labels for nodes are swapped
                atomLabels['*1'].label = '*2'
                atomLabels['*2'].label = '*1'
                atomLabels['*3'].label = '*5'
                atomLabels['*5'].label = '*3'

        if not forward: template = self.reverseTemplate
        else:           template = self.forwardTemplate

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

        # If there are two product structures, place the one containing '*1' first
        if len(productStructures) == 2:
            if not productStructures[0].containsLabeledAtom('*1') and \
                productStructures[1].containsLabeledAtom('*1'):
                productStructures.reverse()

        # If product structures are Molecule objects, update their atom types
        for struct in productStructures:
            if isinstance(struct, Molecule):
                struct.update()

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
        
        productStructures = None

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
        except InvalidActionError:
#            logging.error('Unable to apply reaction recipe!')
#            logging.error('Reaction family is {0} in {1} direction'.format(self.label, 'forward' if forward else 'reverse'))
#            logging.error('Reactant structures are:')
#            for struct in reactantStructures:
#                logging.error(struct.toAdjacencyList())
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

        # If there are two product structures, place the one containing '*1' first
        if len(productStructures) == 2:
            if not productStructures[0].containsLabeledAtom('*1') and \
                productStructures[1].containsLabeledAtom('*1'):
                productStructures.reverse()

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
        from rmgpy.data.rmg import getDB
        
        forbidden_structures = getDB('forbidden')

        if self.forbidden is not None and self.forbidden.isMoleculeForbidden(molecule):
            return True
        if forbidden_structures.isMoleculeForbidden(molecule):
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

        # Convert to species objects
        reactants = [Species(molecule=[mol]) for mol in reactants]
        products = [Species(molecule=[mol]) for mol in products]

        # Create and return template reaction object
        reaction = TemplateReaction(
            reactants = reactants if isForward else products,
            products = products if isForward else reactants,
            degeneracy = 1,
            reversible = True,
            family = self.label,
        )
        
        # Store the labeled atoms so we can recover them later
        # (e.g. for generating reaction pairs and templates)
        labeledAtoms = []
        for reactant in reaction.reactants:
            for label, atom in reactant.molecule[0].getLabeledAtoms().items():
                labeledAtoms.append((label, atom))
        reaction.labeledAtoms = labeledAtoms
        
        return reaction

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

    def generateReactions(self, reactants):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be either single :class:`Molecule` objects
        or lists of same. Does not estimate the kinetics of these reactions
        at this time. Returns a list of :class:`TemplateReaction` objects
        using :class:`Molecule` objects for both reactants and products
        The reactions are constructed such that the forward direction is
        consistent with the template of this reaction family.
        """
        if isinstance(reactants, tuple):
            reactants = list(reactants)
        # We want to work with species objects
        for i, reactant in enumerate(reactants):
            if isinstance(reactant, Molecule):
                reactants[i] = Species(molecule=[reactant])
            elif isinstance(reactant, list):
                reactants[i] = Species(molecule=reactant)
        reactionList = []
        
        # Forward direction (the direction in which kinetics is defined)
        reactionList.extend(self.__generateReactions(reactants, forward=True))
        
        if self.ownReverse:
            # for each reaction, make its reverse reaction and store in a 'reverse' attribute
            for rxn in reactionList:
                reactions = self.__generateReactions(rxn.products, products=rxn.reactants, forward=True)
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
                        reactions = self.__generateReactions(rxn.products, products=rxn.reactants, forward=True)
                    finally:
                        self.forbidden = tempObject
                    if len(reactions) != 1:
                        logging.error("Still experiencing error: Expecting one matching reverse reaction, not {0} in reaction family {1} for forward reaction {2}.\n".format(len(reactions), self.label, str(rxn)))
                        raise KineticsError("Did not find reverse reaction in reaction family {0} for reaction {1}.".format(self.label, str(rxn)))
                    else:
                        logging.error("Error was fixed, the product is a forbidden structure when used as a reactant in the reverse direction.")
                        # Delete this reaction, since it should probably also be forbidden in the initial direction
                        # Hack fix for now
                        del rxn
                elif len(reactions) > 1:
                    logging.error("Expecting one matching reverse reaction, not {0} in reaction family {1} for forward reaction {2}.\n".format(len(reactions), self.label, str(rxn)))
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


            
        else: # family is not ownReverse
            # Reverse direction (the direction in which kinetics is not defined)
            reactionList.extend(self.__generateReactions(reactants, forward=False))
            
        return reactionList
    
    def calculateDegeneracy(self, reaction):
        """
        For a `reaction` given in the direction in which the kinetics are
        defined, compute the reaction-path degeneracy.
        """
        reactions = self.__generateReactions(reaction.reactants, products=reaction.products, forward=True)
        if len(reactions) != 1:
            for reactant in reaction.reactants:
                logging.error(reactant)
            for product in reaction.products:
                logging.error(product)
            raise KineticsError(('Unable to calculate degeneracy for reaction {0} '
                                 'in reaction family {1}. Expected 1 reaction '
                                 'but generated {2}').format(reaction, self.label, len(reactions)))
        return reactions[0].degeneracy
        
    def __generateReactions(self, reactants, products=None, forward=True):
        """
        Generate a list of all of the possible reactions of this family between
        the list of `reactants`. The number of reactants provided must match
        the number of reactants expected by the template, or this function
        will return an empty list. Each item in the list of reactants should
        be a list of :class:`Molecule` objects, each representing a resonance
        isomer of the species of interest.
        """

        for i, reactant in enumerate(reactants):
            if isinstance(reactant, Molecule):
                reactants[i] = Species(molecule=[reactant])
            elif isinstance(reactant, list):
                reactants[i] = Species(molecule=reactant)
        if products is not None:
            for i, product in enumerate(products):
                if isinstance(product, Molecule):
                    products[i] = Species(molecule=[product])
                elif isinstance(product, list):
                    products[i] = Species(molecule=product)

        rxnList = []

        if forward:
            template = self.forwardTemplate
        elif self.reverseTemplate is None:
            return []
        else:
            template = self.reverseTemplate

        # Unimolecular reactants: A --> products
        if len(reactants) == 1 and len(template.reactants) == 1:

            # Iterate over all resonance isomers of the reactant
            for molecule in reactants[0].molecule:

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

            moleculesA = reactants[0].molecule
            moleculesB = reactants[1].molecule

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
                                productStructures = self.__generateProductStructures(reactantStructures, [mapA, mapB], forward)
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
        # If products is given, remove reactions from the reaction list that
        # don't generate the given products
        if products is not None:

            for product in products:
                assert isinstance(product, Species)
                product.generateResonanceIsomers()

            rxnList0 = rxnList[:]
            rxnList = []
            index = 0
            for reaction in rxnList0:
            
                products0 = reaction.products if forward else reaction.reactants

                # For aromatics, generate aromatic resonance structures to accurately identify isomorphic species
                for product0 in products0:
                    product0.molecule.extend(generateAromaticResonanceStructures(product0.molecule[0]))

                # Skip reactions that don't match the given products
                match = False

                if len(products) == len(products0) == 1:
                    if products0[0].isIsomorphic(products[0]):
                        match = True
                elif len(products) == len(products0) == 2:
                    if products0[0].isIsomorphic(products[0]) and products0[1].isIsomorphic(products[1]):
                        match = True
                    elif products0[0].isIsomorphic(products[1]) and products0[1].isIsomorphic(products[0]):
                        match = True

                if match: 
                    rxnList.append(reaction)
            
        # The reaction list may contain duplicates of the same reaction
        # These duplicates should be combined and the reaction degeneracy should be increased
        # Reaction degeneracy is only increased if two reaction are isomorphic but resulted
        # from different transition states.
        # We want to sort all the reactions into sublists composed of isomorphic reactions
        rxnSorted = []
        for rxn0 in rxnList:

            products0 = rxn0.products if forward else rxn0.reactants
            for product in products0:
                assert isinstance(product, Species)
                product.generateResonanceIsomers(keepIsomorphic=True, keepInitial=True, replaceExisting=True)

            if len(rxnSorted) == 0:
                # This is the first reaction, so create a new sublist
                rxnSorted.append([rxn0])
            else:
                # Loop through each sublist, which represents a unique reaction
                for rxnList1 in rxnSorted:
                    # Try to determine if the current rxn0 is identical or isomorphic to any reactions in the sublist
                    identical = False
                    isomorphic = False
                    for rxn in rxnList1:
                        products = rxn.products if forward else rxn.reactants

                        # We know the reactants are the same, so we only need to compare the products
                        if len(products) == len(products0) == 1:
                            if products[0].isIdentical(products0[0]):
                                identical = True
                                isomorphic = True
                            elif not isomorphic and products[0].isIsomorphic(products0[0]):
                                isomorphic = True
                        elif len(products) == len(products0) == 2:
                            if products[0].isIdentical(products0[0]) and products[1].isIdentical(products0[1]):
                                identical = True
                                isomorphic = True
                            elif products[0].isIdentical(products0[1]) and products[1].isIdentical(products0[0]):
                                identical = True
                                isomorphic = True
                            elif not isomorphic and products[0].isIsomorphic(products0[0]) and products[1].isIsomorphic(products0[1]):
                                isomorphic = True
                            elif not isomorphic and products[0].isIsomorphic(products0[1]) and products[1].isIsomorphic(products0[0]):
                                isomorphic = True

                        if identical:
                            # An exact copy of rxn0 is already in our list, so we can move on to the next rxn
                            break
                        elif isomorphic:
                            # This is the right sublist for rxn0, but continue to see if there is an identical rxn
                            continue
                        else:
                            # This is the wrong sublist, so move on to the next sublist
                            break
                    else:
                        # We did not break, so this is the right sublist, but there is no identical reaction
                        # This means that we should add rxn0 to the sublist as a degenerate rxn
                        rxnList1.append(rxn0)
                    if identical or isomorphic:
                        # We already found the right sublist, so we can move on to the next rxn
                        break
                else:
                    # We did not break, which means that there was no isomorphic sublist, so create a new one
                    rxnSorted.append([rxn0])

        rxnList = []
        for rxnList1 in rxnSorted:
            # Collapse our sorted reaction list by taking one reaction from each sublist
            rxn = rxnList1[0]
            # The degeneracy of each reaction is the number of reactions that were in the sublist
            rxn.degeneracy = len(rxnList1)
            rxnList.append(rxn)
        
        # Determine the reactant-product pairs to use for flux analysis
        # Also store the reaction template (useful so we can easily get the kinetics later)
        for reaction in rxnList:
            
            # Restore the labeled atoms long enough to generate some metadata
            for reactant in reaction.reactants:
                for mol in reactant.molecule:
                    mol.clearLabeledAtoms()
            for label, atom in reaction.labeledAtoms:
                atom.label = label
            
            # Generate metadata about the reaction that we will need later
            reaction.pairs = self.getReactionPairs(reaction)
            reaction.template = self.getReactionTemplateLabels(reaction)
            if not forward:
                reaction.degeneracy = self.calculateDegeneracy(reaction)

            # Unlabel the atoms
            for label, atom in reaction.labeledAtoms:
                atom.label = ''
            
            # We're done with the labeled atoms, so delete the attribute
            del reaction.labeledAtoms
            
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
            if reaction.reactants[0].molecule[0].containsLabeledAtom('*1'):
                if reaction.products[0].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
            elif reaction.reactants[1].molecule[0].containsLabeledAtom('*1'):
                if reaction.products[1].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
        elif self.label.lower() == 'disproportionation':
            # Hardcoding for disproportionation: pair the reactant containing
            # *1 with the product containing *1
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].molecule[0].containsLabeledAtom('*1'):
                if reaction.products[0].molecule[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].molecule[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
            elif reaction.reactants[1].molecule[0].containsLabeledAtom('*1'):
                if reaction.products[1].molecule[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].molecule[0].containsLabeledAtom('*1'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
        elif self.label.lower() in ['substitution_o', 'substitutions']:
            # Hardcoding for Substitution_O: pair the reactant containing
            # *2 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].molecule[0].containsLabeledAtom('*2'):
                if reaction.products[0].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[1].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[1]])
                    pairs.append([reaction.reactants[1],reaction.products[0]])
                else:
                    error = True
            elif reaction.reactants[1].molecule[0].containsLabeledAtom('*2'):
                if reaction.products[1].molecule[0].containsLabeledAtom('*3'):
                    pairs.append([reaction.reactants[0],reaction.products[0]])
                    pairs.append([reaction.reactants[1],reaction.products[1]])
                elif reaction.products[0].molecule[0].containsLabeledAtom('*3'):
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
        for kinetics, entry, isForward in kineticsList:
            if kinetics is not None:
                kinetics.comment += "Matched reaction {0} {1} in {2}".format(entry.index, entry.label, depository.label)
        return kineticsList
    
    def __selectBestKinetics(self, kineticsList):
        """
        For a given set of kinetics `kineticsList`, return the kinetics deemed
        to be the "best". This is determined to be the one with the lowest
        non-zero rank that occurs first.
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
        results, with each result comprising the kinetics, the source, and
        the entry. If it came from a template estimate, the source and entry
        will both be `None`.
        If returnAllKinetics==False, only the first (best?) matching kinetics is returned.
        """
        kineticsList = []
        
        depositories = self.depositories[:]

        template = self.retrieveTemplate(templateLabels)
        
        # Check the various depositories for kinetics
        for depository in depositories:
            kineticsList0 = self.getKineticsFromDepository(depository, reaction, template, degeneracy)
            if len(kineticsList0) > 0 and not returnAllKinetics:
                kinetics, entry, isForward = self.__selectBestKinetics(kineticsList0)
                return kinetics, depository, entry, isForward
            else:
                for kinetics, entry, isForward in kineticsList0:
                    kineticsList.append([kinetics, depository, entry, isForward])
        
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
        
    def estimateKineticsUsingRateRules(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.
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
            for map in mappings:
                reactantStructures = [molecule]
                try:
                    productStructures = self.__generateProductStructures(reactantStructures, [map], forward=True)
                except ForbiddenStructureException:
                    pass
                else:
                    if productStructures is not None:
                        if len(products) == 1 and len(productStructures) == 1:
                            if products[0].isIsomorphic(productStructures[0]):
                                return reactantStructures, productStructures
                        elif len(products) == 2 and len(productStructures) == 2:
                            if products[0].isIsomorphic(productStructures[0]):
                                if products[1].isIsomorphic(productStructures[1]):
                                    return reactantStructures, productStructures
                            if products[0].isIsomorphic(productStructures[1]):
                                if products[1].isIsomorphic(productStructures[0]):
                                    return reactantStructures, productStructures
                        else: continue
            # if there're some mapping available but cannot match the provided products
            # raise exception
            if len(mappings) > 0:
                raise Exception('Something wrong with products that RMG cannot find a match!')
            return (None, None)
        elif len(reactants0) == 2:
            moleculeA = reactants0[0]
            moleculeB = reactants0[1]
            mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[0])
            mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[1])

            # Iterate over each pair of matches (A, B)
            for mapA in mappingsA:
                for mapB in mappingsB:
                    reactantStructures = [moleculeA, moleculeB]
                    try:
                        productStructures = self.__generateProductStructures(reactantStructures, [mapA, mapB], forward=True)
                    except ForbiddenStructureException:
                        pass
                    else:
                        if productStructures is not None:
                            if len(products) == 1 and len(productStructures) == 1:
                                if products[0].isIsomorphic(productStructures[0]):
                                    return reactantStructures, productStructures
                            elif len(products) == 2 and len(productStructures) == 2:
                                if products[0].isIsomorphic(productStructures[0]):
                                    if products[1].isIsomorphic(productStructures[1]):
                                        return reactantStructures, productStructures
                                if products[0].isIsomorphic(productStructures[1]):
                                    if products[1].isIsomorphic(productStructures[0]):
                                        return reactantStructures, productStructures
                            else: continue
            # if there're some mapping available but cannot match the provided products
            # raise exception
            if len(mappingsA)*len(mappingsB) > 0:
                raise Exception('Something wrong with products that RMG cannot find a match!')
            return (None, None)
        else:
            raise Exception('You have {0} reactants, which is unexpected!'.format(len(reactants)))
        
    def addAtomLabelsForReaction(self, reaction):
        """
        Apply atom labels on a reaction using the appropriate atom labels from this this reaction family.  
        The reaction's reactants and products must be lists of Molecule objects.
        The reaction is modified to use Species objects containing the labeled molecules after this function.
        """
        reactants = reaction.reactants[:]
        products = reaction.products[:]
       
        labeledReactants, labeledProducts = self.getLabeledReactantsAndProducts(reactants, products)
        if not labeledReactants and not labeledProducts:
            if len(reactants) > 1:
                # Check for swapped reactants if there are more than 1
                labeledReactants, labeledProducts = self.getLabeledReactantsAndProducts(list(reversed(reactants)), products)
        if not labeledReactants and not labeledProducts:
            raise Exception("RMG could not label atoms for this reaction in {}.".format(self.label))
        
        labeledProducts_spcs = []
        labeledReactants_spcs = []
        for labeledProduct in labeledProducts:
            labeledProducts_spcs.append(Species(molecule=[labeledProduct]))
        reaction.products = labeledProducts_spcs
        for labeledReactant in labeledReactants:
            labeledReactants_spcs.append(Species(molecule=[labeledReactant]))
        reaction.reactants = labeledReactants_spcs

    def getTrainingDepository(self):
        """
        Returns the `training` depository from self.depositories
        """
        for depository in self.depositories:
            if depository.label.endswith('training'):
                return depository
        else:
            raise Exception('Could not find training depository in family {0}.'.format(self.label))
        
    def retrieveOriginalEntry(self, templateLabel):
        """
        Retrieves the original entry, be it a rule or training reaction, given
        the template label in the form 'group1;group2' or 'group1;group2;group3'
        
        Returns tuple in the form
        (RateRuleEntry, TrainingReactionEntry)
        
        Where the TrainingReactionEntry is only present if it comes from a training reaction
        """
        
        templateLabels = templateLabel.split()[0].split(';')
        template = self.retrieveTemplate(templateLabels)
        rule = self.getRateRule(template)
        if 'from training reaction' in rule.data.comment:
            trainingIndex = int(rule.data.comment.split()[-1])
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
            # Discard the last line, unless it's the only line!
            # The last line is 'Estimated using ... for rate rule (originalTemplate)'
            if len(lines) == 1:
                comment = lines[0]
                if comment.startswith('Estimated using template'):
                    tokenTemplateLabel = comment.split()[3][1:-1]
                    ruleEntry, trainingEntry = self.retrieveOriginalEntry(tokenTemplateLabel) 
                    if trainingEntry:
                        training.append((ruleEntry,trainingEntry,1))   # Weight is 1
                    else:
                        rules.append((ruleEntry,1))
                else:
                    raise Exception('Could not parse unexpected line found in kinetics comment: {}'.format(comment))
            else:
                comment = ' '.join(lines[:-1])
                # Clean up line for exec
                evalCommentString = re.sub(r" \+ ", ",",                        # any remaining + signs
                                    re.sub(r"Average of ", "",                  # average of averages
                                    re.sub(r"Average of \[(?!Average)", "['",   # average of groups
                                    re.sub(r"(\b|\))]", r"\1']",                # initial closing bracket
                                    re.sub(r"(?<=\b) \+ (?=Average)", "',",     # + sign between non-average and average
                                    re.sub(r"(?<=]) \+ (?!Average)", ",'",      # + sign between average and non-average
                                    re.sub(r"(?<!]) \+ (?!Average)", "','",     # + sign between non-averages
                                    comment)))))))

                entryNestedList = eval(evalCommentString)
                
                weightedEntries = assignWeightsToEntries(entryNestedList, [])
                
                
                rules = {}
                training = {}
                
                for tokenTemplateLabel, weight in weightedEntries:
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
        
        For a training reaction, the Source_Data returns 
        [Family_Label, Training_Reaction_Entry, Kinetics_In_Reverse?]
        
        For a reaction from rate rules, the Source_Data is a tuple containing
        [Family_Label, {'template': originalTemplate,
            'degeneracy': degeneracy, 
            'exact': boolean_exact?, 
            'rules': a list of (original rate rule entry, weight in average)
            'training': a list of (original rate rule entry associated with training entry, original training entry, weight in average)
            }]
        where TrainingReactions are ones that have created rules used in the estimate.
        
        where Exact is a boolean of whether the rate is an exact match, 
        Template is the reaction template used,
        and RateRules is a list of the rate rule entries containing the kinetics used
        """
        import re
        lines = reaction.kinetics.comment.split('\n')

        exact = False
        template = None
        rules = None
        trainingEntries = None
        degeneracy = 1

        regex = "\((.*)\)" # only hit outermost parentheses
        for line in lines:
            if line.startswith('Matched'):
                # Source of the kinetics is from training reaction
                trainingReactionIndex = int(line.split()[2])
                depository  = self.getTrainingDepository()
                trainingEntry = depository.entries[trainingReactionIndex]
                # Perform sanity check that the training reaction's label matches that of the comments
                if trainingEntry.label not in line:
                    raise Exception('Reaction {0} uses kinetics from training reaction {1} but does not match the training reaction {1} from the {2} family.'.format(reaction,trainingReactionIndex,self.label))
                
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
                degeneracy = int(line.split()[-1])

        # Extract the rate rule information 
        fullCommentString = reaction.kinetics.comment.replace('\n', ' ')
        
        # The rate rule string is right after the phrase 'for rate rule'
        rateRuleString = fullCommentString.split("for rate rule",1)[1].split()[0]
        templateLabel = re.split(regex, rateRuleString)[1]
        template = self.retrieveTemplate(templateLabel.split(';'))
        rules, trainingEntries = self.getSourcesForTemplate(template)
        

        if not template:
            raise Exception('Could not extract kinetics source from comments for reaction {}.'.format(reaction))
        
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

