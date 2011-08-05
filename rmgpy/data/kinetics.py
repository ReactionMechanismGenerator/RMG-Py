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
import re
import logging
import codecs
from copy import copy, deepcopy

from base import *

from rmgpy.quantity import Quantity
from rmgpy.reaction import Reaction, ReactionError
from rmgpy.kinetics import *
from rmgpy.group import GroupBond, Group
from rmgpy.molecule import Bond
from rmgpy.species import Species

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

class DepositoryReaction(Reaction):
    """
    A Reaction object generated from a reaction depository. In addition to the
    usual attributes, this class includes `depository` and `entry` attributes to
    store the library and the entry in that depository that it was created from.
    """

    def __init__(self, index=-1, reactants=None, products=None, kinetics=None, reversible=True, transitionState=None, thirdBody=False, duplicate=False, degeneracy=1, depository=None, family=None, entry=None):
        Reaction.__init__(self, index=index, reactants=reactants, products=products, kinetics=kinetics, reversible=reversible, transitionState=transitionState, thirdBody=thirdBody, duplicate=duplicate, degeneracy=degeneracy)
        self.depository = depository
        self.family = family
        self.entry = entry

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (DepositoryReaction, (self.index, self.reactants, self.products, self.kinetics, self.reversible, self.transitionState, self.thirdBody, self.degeneracy, self.depository, self.entry))

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        DepositoryReaction this should be a KineticsDepository object.
        """
        return self.depository

################################################################################

class LibraryReaction(Reaction):
    """
    A Reaction object generated from a reaction library. In addition to the
    usual attributes, this class includes `library` and `entry` attributes to
    store the library and the entry in that library that it was created from.
    """
    
    def __init__(self, index=-1, reactants=None, products=None, kinetics=None, reversible=True, transitionState=None, thirdBody=False, duplicate=False, degeneracy=1, library=None, entry=None):
        Reaction.__init__(self, index=index, reactants=reactants, products=products, kinetics=kinetics, reversible=reversible, transitionState=transitionState, thirdBody=thirdBody, duplicate=duplicate, degeneracy=degeneracy)
        self.library = library
        self.family = library
        self.entry = entry

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (LibraryReaction, (self.index, self.reactants, self.products, self.kinetics, self.reversible, self.transitionState, self.thirdBody, self.degeneracy, self.library, self.entry))

    def getSource(self):
        """
        Return the database that was the source of this reaction. For a
        LibraryReaction this should be a KineticsLibrary object.
        """
        return self.library

################################################################################

class TemplateReaction(Reaction):
    """
    A Reaction object generated from a reaction family template. In addition to
    the usual attributes, this class includes a `family` attribute to store the
    family that it was created from.
    """

    def __init__(self, index=-1, reactants=None, products=None, kinetics=None, reversible=True, transitionState=None, thirdBody=False, duplicate=False, degeneracy=1, family=None, template=None):
        Reaction.__init__(self, index=index, reactants=reactants, products=products, kinetics=kinetics, reversible=reversible, transitionState=transitionState, thirdBody=thirdBody, duplicate=duplicate, degeneracy=degeneracy)
        self.family = family
        self.template = template

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TemplateReaction, (self.index, self.reactants, self.products, self.kinetics, self.reversible, self.transitionState, self.thirdBody, self.degeneracy, self.family, self.template))

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
                    bond = GroupBond(order=['S']) if pattern else Bond(order='S')
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

def saveEntry(f, entry):
    """
    Save an `entry` in the kinetics database by writing a string to
    the given file object `f`.
    """

    def sortEfficiencies(efficiencies0):
        efficiencies = {}
        for mol, eff in efficiencies0.iteritems():
            smiles = mol.toSMILES()
            efficiencies[smiles] = eff
        keys = efficiencies.keys()
        keys.sort()
        return [(key, efficiencies[key]) for key in keys]

    f.write('entry(\n')
    f.write('    index = {0:d},\n'.format(entry.index))
    if entry.label != '':
        f.write('    label = "{0}",\n'.format(entry.label))

    if isinstance(entry.item, Reaction):
        for i, reactant in enumerate(entry.item.reactants):
            if isinstance(reactant, Molecule):
                f.write('    reactant{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(reactant.toAdjacencyList(removeH=True))
                f.write('""",\n')
            elif isinstance(reactant, Species):
                f.write('    reactant{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(reactant.molecule[0].toAdjacencyList(label=reactant.label, removeH=True))
                f.write('""",\n')
            elif isinstance(reactant, Group):
                f.write('    group{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(reactant.toAdjacencyList())
                f.write('""",\n')
            elif isinstance(reactant, LogicNode):
                f.write('    group{0:d} = "{1}",\n'.format(i+1, reactant))
        for i, product in enumerate(entry.item.products):
            if isinstance(product, Molecule):
                f.write('    product{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(product.toAdjacencyList(removeH=True))
                f.write('""",\n')
            elif isinstance(reactant, Species):
                f.write('    product{0:d} = \n'.format(i+1))
                f.write('"""\n')
                f.write(product.molecule[0].toAdjacencyList(label=product.label, removeH=True))
                f.write('""",\n')
        if not isinstance(entry.item.reactants[0], Group) and not isinstance(entry.item.reactants[0], LogicNode):
            f.write('    degeneracy = {0:d},\n'.format(entry.item.degeneracy))
        if entry.item.duplicate: 
            f.write('    duplicate = {0!r},\n'.format(entry.item.duplicate))
    elif isinstance(entry.item, Group):
        f.write('    group = \n')
        f.write('"""\n')
        f.write(entry.item.toAdjacencyList())
        f.write('""",\n')
    elif isinstance(entry.item, LogicNode):
        f.write('    group = "{0}",\n'.format(entry.item))
    else:
        raise InvalidDatabaseError("Encountered unexpected item of type {0} while saving database.".format(entry.item.__class__))

    if isinstance(entry.data, Arrhenius):
        f.write('    kinetics = Arrhenius(\n')
        f.write('        A = {0!r},\n'.format(entry.data.A))
        f.write('        n = {0!r},\n'.format(entry.data.n))
        f.write('        Ea = {0!r},\n'.format(entry.data.Ea))
        f.write('        T0 = {0!r},\n'.format(entry.data.T0))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, ArrheniusEP):
        f.write('    kinetics = ArrheniusEP(\n')
        f.write('        A = {0!r},\n'.format(entry.data.A))
        f.write('        n = {0!r},\n'.format(entry.data.n))
        f.write('        alpha = {0!r},\n'.format (entry.data.alpha))
        f.write('        E0 = {0!r},\n'.format(entry.data.E0))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        f.write('    ),\n')
    elif isinstance(entry.data, MultiKinetics):
        f.write('    kinetics = MultiKinetics(\n')
        f.write('        kineticsList = [')
        for kin in entry.data.kineticsList:
            for line in repr(kin).split('\n'):
                f.write('\n            '+line)
            if kin is entry.data.kineticsList[-1]:
                f.write('\n')
            else:
                f.write(',')
        f.write('        ],\n')
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.Pmin is not None: f.write('        Pmin = {0!r},\n'.format(entry.data.Pmin))
        if entry.data.Pmax is not None: f.write('        Pmax = {0!r},\n'.format(entry.data.Pmax))
        f.write('    ),\n')
    elif isinstance(entry.data, PDepArrhenius):
        # Why is all this here AND in rmgpy.kinetics.PDepArrhenius.__repr__ ?
        f.write('    kinetics = PDepArrhenius(\n')
        f.write('        pressures = {0!r},\n'.format(entry.data.pressures))
        f.write('        arrhenius = [\n')
        for arrh in entry.data.arrhenius:
            f.write('            {0!r},\n'.format(arrh))
        f.write('        ],\n')
        if entry.data.highPlimit is not None: f.write('        highPlimit = {0!r},\n'.format(entry.data.highPlimit))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.Pmin is not None: f.write('        Pmin = {0!r},\n'.format(entry.data.Pmin))
        if entry.data.Pmax is not None: f.write('        Pmax = {0!r},\n'.format(entry.data.Pmax))
        f.write('    ),\n')
    elif isinstance(entry.data, Troe):
        f.write('    kinetics = Troe(\n')
        f.write('        arrheniusHigh = {0!r},\n'.format(entry.data.arrheniusHigh))
        f.write('        arrheniusLow = {0!r},\n'.format(entry.data.arrheniusLow))
        f.write('        efficiencies = {{{0}}},\n'.format(', '.join(['"{0}": {1:g}'.format(label, eff) for label, eff in sortEfficiencies(entry.data.efficiencies)])))
        f.write('        alpha = {0!r},\n'.format(entry.data.alpha))
        f.write('        T3 = {0!r},\n'.format(entry.data.T3))
        f.write('        T1 = {0!r},\n'.format(entry.data.T1))
        if entry.data.T2 is not None: f.write('        T2 = {0!r},\n'.format(entry.data.T2))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.Pmin is not None: f.write('        Pmin = {0!r},\n'.format(entry.data.Pmin))
        if entry.data.Pmax is not None: f.write('        Pmax = {0!r},\n'.format(entry.data.Pmax))
        f.write('    ),\n')
    elif isinstance(entry.data, Lindemann):
        f.write('    kinetics = Lindemann(\n')
        f.write('        arrheniusHigh = {0!r},\n'.format(entry.data.arrheniusHigh))
        f.write('        arrheniusLow = {0!r},\n'.format(entry.data.arrheniusLow))
        f.write('        efficiencies = {{{0}}},\n'.format(', '.join(['"{0}": {1:g}'.format(label, eff) for label, eff in sortEfficiencies(entry.data.efficiencies)])))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.Pmin is not None: f.write('        Pmin = {0!r},\n'.format(entry.data.Pmin))
        if entry.data.Pmax is not None: f.write('        Pmax = {0!r},\n'.format(entry.data.Pmax))
        f.write('    ),\n')
    elif isinstance(entry.data, ThirdBody):
        f.write('    kinetics = ThirdBody(\n')
        f.write('        arrheniusHigh = {0!r},\n'.format(entry.data.arrheniusHigh))
        f.write('        efficiencies = {{{0}}},\n'.format(', '.join(['"{0}": {1:g}'.format(label, eff) for label, eff in sortEfficiencies(entry.data.efficiencies)])))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.Pmin is not None: f.write('        Pmin = {0!r},\n'.format(entry.data.Pmin))
        if entry.data.Pmax is not None: f.write('        Pmax = {0!r},\n'.format(entry.data.Pmax))
        f.write('    ),\n')
    elif isinstance(entry.data, Chebyshev):
        f.write('    kinetics = Chebyshev(\n')
        f.write('        coeffs = [\n')
        for i in range(entry.data.degreeT):
            f.write('            [{0}],\n'.format(','.join(['{0:g}'.format(self.coeffs[i,j]) for j in range(self.degreeP)])))
        f.write('        ],\n')
        f.write('        kunits = {0},\n'.format(entry.data.kunits))
        if entry.data.Tmin is not None: f.write('        Tmin = {0!r},\n'.format(entry.data.Tmin))
        if entry.data.Tmax is not None: f.write('        Tmax = {0!r},\n'.format(entry.data.Tmax))
        if entry.data.Pmin is not None: f.write('        Pmin = {0!r},\n'.format(entry.data.Pmin))
        if entry.data.Pmax is not None: f.write('        Pmax = {0!r},\n'.format(entry.data.Pmax))
        f.write('    ),\n')
    else:
        f.write('    kinetics = {0!r},\n'.format(entry.data))

    f.write('    reference = {0!r},\n'.format(entry.reference))
    f.write('    referenceType = "{0}",\n'.format(entry.referenceType))
    if entry.rank is not None:
        f.write('    rank = {0},\n'.format(entry.rank))
    f.write('    shortDesc = u"""{0}""",\n'.format(entry.shortDesc))
    f.write('    longDesc = \n')
    f.write('u"""\n')
    f.write(entry.longDesc.strip() + "\n")
    f.write('""",\n')

    f.write('    history = [\n')
    for time, user, action, description in entry.history:
        f.write('        ("{0}","{1}","{2}","""{3}"""),\n'.format(time, user, action, description))
    f.write('    ],\n')

    f.write(')\n\n')

################################################################################

class KineticsDepository(Database):
    """
    A class for working with an RMG kinetics depository. Each depository 
    corresponds to a reaction family (a :class:`KineticsGroups` object). Each
    entry in a kinetics depository involves a reaction defined either by
    real reactant and product species (as in a kinetics library) or a set of
    functional groups (as in a reaction family).
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def __repr__(self):
        return '<KineticsDepository "{0}">'.format(self.label)

    def loadEntry(self, index, reactant1=None, reactant2=None, reactant3=None, product1=None, product2=None, product3=None, group1=None, group2=None, group3=None, kinetics=None, degeneracy=1, label='', reference=None, referenceType='', shortDesc='', longDesc='', rank=None, history=None):
        
        if reactant1 is not None and product1 is not None:
            # The reaction involves real reactants and products
            assert group1 is None and group2 is None and group3 is None
            
            reactants = [Molecule().fromAdjacencyList(reactant1)]
            if reactant2 is not None: reactants.append(Molecule().fromAdjacencyList(reactant2))
            if reactant3 is not None: reactants.append(Molecule().fromAdjacencyList(reactant3))

            products = [Molecule().fromAdjacencyList(product1)]
            if product2 is not None: products.append(Molecule().fromAdjacencyList(product2))
            if product3 is not None: products.append(Molecule().fromAdjacencyList(product3))
            
            reaction = Reaction(reactants=reactants, products=products, degeneracy=degeneracy)
        
        elif group1 is not None:
            # The reaction involves functional groups
            assert reactant1 is None and reactant2 is None and reactant3 is None
            assert product1 is None and product2 is None and product3 is None
            
            reactants = []
            
            if group1[0:3].upper() == 'OR{' or group1[0:4].upper() == 'AND{' or group1[0:7].upper() == 'NOT OR{' or group1[0:8].upper() == 'NOT AND{':
                reactants.append(makeLogicNode(group1))
            else:
                reactants.append(Group().fromAdjacencyList(group1))
            if group2 is not None: 
                if group2[0:3].upper() == 'OR{' or group2[0:4].upper() == 'AND{' or group2[0:7].upper() == 'NOT OR{' or group2[0:8].upper() == 'NOT AND{':
                    reactants.append(makeLogicNode(group2))
                else:
                    reactants.append(Group().fromAdjacencyList(group2))
            if group3 is not None: 
                if group3[0:3].upper() == 'OR{' or group3[0:4].upper() == 'AND{' or group3[0:7].upper() == 'NOT OR{' or group3[0:8].upper() == 'NOT AND{':
                    reactants.append(makeLogicNode(group3))
                else:
                    reactants.append(Group().fromAdjacencyList(group3))
            
            reaction = Reaction(reactants=reactants, products=[])
        else:
            raise ValueError("Entry doesn't seem to involve reactants and products, or groups.")
            
        entry = Entry(
            index = index,
            label = label,
            item = reaction,
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            rank = rank,
            history = history or [],
        )
        self.entries[index] = entry
        return entry

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return saveEntry(f, entry)

    def processOldLibraryEntry(self, data):
        """
        Process a list of parameters `data` as read from an old-style RMG
        thermo database, returning the corresponding kinetics object.
        """
        # This is hardcoding of reaction families!
        label = os.path.split(self.label)[-2]
        if label in ['H_Abstraction', 'R_Addition_MultipleBond', 'R_Recombination', 'HO2_Elimination_from_PeroxyRadical', 'Disproportionation', '1+2_Cycloaddition', '2+2_cycloaddition_Cd', '2+2_cycloaddition_CO', '2+2_cycloaddition_CCO', 'Diels_alder_addition', '1,2_Insertion', '1,3_Insertion_CO2', '1,3_Insertion_ROR', 'R_Addition_COm', 'Oa_R_Recombination']:
            Aunits = 'cm^3/(mol*s)'
        elif label in ['intra_H_migration', 'Birad_recombination', 'intra_OH_migration', 'Cyclic_Ether_Formation', 'Intra_R_Add_Exocyclic', 'Intra_R_Add_Endocyclic', '1,2-Birad_to_alkene', 'Intra_Disproportionation']:
            Aunits = 's^-1'
        else:
            raise ValueError('Unable to determine preexponential units for old reaction family "{0}".'.format(self.label))

        try:
            Tmin, Tmax = data[0].split('-')
            Tmin = (float(Tmin),"K")
            Tmax = (float(Tmax),"K")
        except ValueError:
            Tmin = (float(data[0]),"K")
            Tmax = None

        A, n, alpha, E0, dA, dn, dalpha, dE0 = data[1:9]
        
        A = float(A)
        if dA[0] == '*':
            A = Quantity(A,Aunits,'*|/',float(dA[1:]))
        else:
            dA = float(dA)
            if dA != 0:
                A = Quantity(A,Aunits,'+|-',dA)
            else:
                A = Quantity(A,Aunits)
        
        n = float(n); dn = float(dn)
        if dn != 0:
            n = Quantity(n,'','+|-',dn)
        else:
            n = Quantity(n,'')
                
        alpha = float(alpha); dalpha = float(dalpha)
        if dalpha != 0:
            alpha = Quantity(alpha,'','+|-',dalpha)
        else:
            alpha = Quantity(alpha,'')
        
        E0 = float(E0); dE0 = float(dE0)
        if dE0 != 0:
            E0 = Quantity(E0,'kcal/mol','+|-',dE0)
        else:
            E0 = Quantity(E0,'kcal/mol')
        
        rank = int(data[9])
        
        return ArrheniusEP(A=A, n=n, alpha=alpha, E0=E0, Tmin=Tmin, Tmax=Tmax), rank

    def loadOldRateRules(self, path, groups):
        """
        Load a set of old rate rules for kinetics groups into this depository.
        """
        # Parse the old library
        numLabels = max(len(groups.forwardTemplate.reactants), len(groups.top))
        entries = self.parseOldLibrary(os.path.join(path, 'rateLibrary.txt'), numParameters=10, numLabels=numLabels)
        
        self.entries = {}
        for entry in entries:
            index, label, data, shortDesc = entry
            kinetics, rank = data
            reactants = [groups.entries[l].item for l in label.split(';')]
            item = Reaction(reactants=reactants, products=[])
            self.entries[index] = Entry(
                index = index,
                label = label,
                item = item,
                data = kinetics,
                rank = rank,
                shortDesc = shortDesc
            )
        self.__loadOldComments(path)
    
    def __loadOldComments(self, path):
        """
        Load a set of old comments from the ``comments.txt`` file for the old
        kinetics groups. This function assumes that the groups have already
        been loaded.
        """
        index = 'General' #mops up comments before the first rate ID
        
        re_underline = re.compile('^\-+')
        
        comments = {}
        comments[index] = ''
        
        # Load the comments into a temporary dictionary for now
        # If no comments file then do nothing
        try:
            f = codecs.open(os.path.join(path, 'comments.rst'), 'r', 'utf-8')
        except IOError:
            return
        for line in f:
            match = re_underline.match(line)
            if match:
                index = f.next().strip()
                assert line.rstrip() == f.next().rstrip(), "Overline didn't match underline"
                if not comments.has_key(index):
                    comments[index] = ''
                line = f.next()
            comments[index] += line
        f.close()
        
        # Transfer the comments to the longDesc attribute of the associated entry
        unused = []
        for index, longDesc in comments.iteritems():
            try:
                index = int(index)
            except ValueError:
                unused.append(index)
                
            if isinstance(index, int):
                for entry in self.entries.values():
                    if entry.index == index:
                        entry.longDesc = longDesc
                        break
                #else:
                #    unused.append(str(index))
            
        # Any unused comments are placed in the longDesc attribute of the depository
        self.longDesc = comments['General'] + '\n'
        unused.remove('General')
        for index in unused:
            try:
                self.longDesc += comments[index] + '\n'
            except KeyError:
                import pdb; pdb.set_trace()
                
    def saveOldRateRules(self, path, groups):
        """
        Save a set of old rate rules for kinetics groups from this depository.
        """
        
        # This is hardcoding of reaction families!
        label = os.path.split(self.label)[-2]
        if label in ['H_Abstraction', 'R_Addition_MultipleBond', 'R_Recombination', 'HO2_Elimination_from_PeroxyRadical', 'Disproportionation', '1+2_Cycloaddition', '2+2_cycloaddition_Cd', '2+2_cycloaddition_CO', '2+2_cycloaddition_CCO', 'Diels_alder_addition', '1,2_Insertion', '1,3_Insertion_CO2', '1,3_Insertion_ROR', 'R_Addition_COm', 'Oa_R_Recombination']:
            factor = 1.0e6
        elif label in ['intra_H_migration', 'Birad_recombination', 'intra_OH_migration', 'Cyclic_Ether_Formation', 'Intra_R_Add_Exocyclic', 'Intra_R_Add_Endocyclic', '1,2-Birad_to_alkene', 'Intra_Disproportionation']:
            factor = 1.0
        else:
            raise ValueError('Unable to determine preexponential units for old reaction family "{0}".'.format(self.label))

        entries = self.entries.values()
        entries.sort(key=lambda x: x.index)
        
        flib = open(os.path.join(path, 'rateLibrary.txt'), 'w')
        flib.write('// The format for the data in this rate library\n')
        flib.write('Arrhenius_EP\n\n')
        
        fcom = codecs.open(os.path.join(path, 'comments.rst'), 'w', 'utf-8')
        fcom.write('-------\n')
        fcom.write('General\n')
        fcom.write('-------\n')
        fcom.write(self.longDesc.strip() + '\n\n')
        
        for entry in entries:
            flib.write('{0:<5d} '.format(entry.index))
            for label in entry.label.split(';'):
                flib.write('{0:<23} '.format(label))
            if entry.data.Tmax is None:
                Trange = '{0:g}    '.format(entry.data.Tmin.value)
            else:
                Trange = '{0:g}-{1:g}    '.format(entry.data.Tmin.value, entry.data.Tmax.value)
            flib.write('{0:<12}'.format(Trange))
            flib.write('{0:11.2e} {1:9.2f} {2:9.2f} {3:11.2f} '.format(entry.data.A.value * factor, entry.data.n.value, entry.data.alpha.value, entry.data.E0.value / 4184.))
            if entry.data.A.isUncertaintyMultiplicative():
                flib.write('*{0:<6g} '.format(entry.data.A.uncertainty))
            else:
                flib.write('{0:<7g} '.format(entry.data.A.uncertainty * factor))
            flib.write('{0:6g} {1:6g} {2:6g} '.format(entry.data.n.uncertainty, entry.data.alpha.uncertainty, entry.data.E0.uncertainty / 4184.))
            flib.write('    {0:<4d}     {1}\n'.format(entry.rank, entry.shortDesc))
            
            fcom.write('------\n')
            fcom.write('{0}\n'.format(entry.index))
            fcom.write('------\n')
            fcom.write(entry.longDesc.strip() + '\n\n')
            
        flib.close()
        fcom.close()
    
################################################################################

class KineticsLibrary(Database):
    """
    A class for working with an RMG kinetics library.
    """

    def __init__(self, label='', name='', shortDesc='', longDesc=''):
        Database.__init__(self, label=label, name=name, shortDesc=shortDesc, longDesc=longDesc)

    def __repr__(self):
        return '<KineticsLibrary "{0}">'.format(self.label)

    def getSpecies(self):
        """
        Return a dictionary containing all of the species in this kinetics
        library.
        """
        speciesDict = {}
        
        def speciesMatch(speciesA, speciesB):
            for moleculeA in speciesA.molecule:
                for moleculeB in speciesB.molecule:
                    if moleculeA.isIsomorphic(moleculeB):
                        return True
            return False
        
        entries = self.entries.values()
        for entry in entries:
            for reactant in entry.item.reactants:
                if reactant.label not in speciesDict:
                    speciesDict[reactant.label] = reactant
                elif not speciesMatch(reactant, speciesDict[reactant.label]):
                    print reactant.molecule[0].toAdjacencyList()
                    print speciesDict[reactant.label].molecule[0].toAdjacencyList()
                    raise DatabaseError('Species label "{0}" used for multiple species in kinetics library {1}.'.format(reactant.label, self.label))
            for product in entry.item.products:
                if product.label not in speciesDict:
                    speciesDict[product.label] = product
                elif not speciesMatch(product, speciesDict[product.label]):
                    import pdb; pdb.set_trace()
                    print product.molecule[0].toAdjacencyList()
                    print speciesDict[product.label].molecule[0].toAdjacencyList()
                    print product.molecule[0].isIsomorphic(speciesDict[product.label].molecule[0])
                    raise DatabaseError('Species label "{0}" used for multiple species in kinetics library {1}.'.format(product.label, self.label))
        
        return speciesDict
    
    def checkForDuplicates(self):
        """
        Check that all duplicate reactions in the kinetics library are
        properly marked (i.e. with their ``duplicate`` attribute set to 
        ``True``).
        """
        for entry0 in self.entries.values():
            reaction0 = entry0.item
            if not reaction0.duplicate:
                # This reaction is not marked as a duplicate reaction
                # This means that if we find any duplicate reactions, it is an error
                for entry in self.entries.values():
                    reaction = entry.item
                    if reaction0 is not reaction and reaction0.reactants == reaction.reactants and reaction0.products == reaction.products:
                        # We found a duplicate reaction that wasn't marked!
                        raise DatabaseError('Unexpected duplicate reaction {0} in kinetics library {1}.'.format(reaction0, self.label))                   

    def convertDuplicatesToMulti(self):
        """
        Merge all marked duplicate reactions in the kinetics library
        into single reactions with multiple kinetics.
        """
        print "trying to find duplicates"
        entries_to_remove = []
        for entry0 in self.entries.values():
            if entry0 in entries_to_remove:
                continue
            reaction0 = entry0.item
            if not reaction0.duplicate:
                continue
            print "Found a duplicate reaction: {0}".format(reaction0)
            duplicates = [entry0]
            for entry in self.entries.values():
                reaction = entry.item
                if reaction0 is reaction:
                    continue
                if reaction0.isIsomorphic(reaction, eitherDirection=False):
                    duplicates.append(entry)
            
            assert len(duplicates)>1
            kineticsList = []
            longDesc = ''
            
            for entry in duplicates:
                kinetics = entry.data
                kineticsList.append(kinetics)
                Tmin = kinetics.Tmin
                Tmax = kinetics.Tmax
                Pmin = kinetics.Pmin
                Pmax = kinetics.Pmax
                longDesc += entry.longDesc+'\n'
            
            entry0.data = MultiKinetics(kineticsList=kineticsList, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax)
            entry0.longDesc = longDesc
            entries_to_remove.extend(duplicates[1:])
        for entry in entries_to_remove:
            print "removing duplicate reaction with index {0}.".format(entry.index)
            del(self.entries[entry.index])
        print "NB. the entries have not been renumbered, so these indices are missing."
        
        
    def load(self, path, local_context=None, global_context=None):
        Database.load(self, path, local_context, global_context)
        
        # Generate a unique set of the species in the kinetics library
        speciesDict = self.getSpecies()
        # Make sure all of the reactions draw from only this set
        entries = self.entries.values()
        for entry in entries:
            entry.item.reactants = [speciesDict[spec.label] for spec in entry.item.reactants]
            entry.item.products = [speciesDict[spec.label] for spec in entry.item.products]
            
        self.checkForDuplicates()
        
    def loadEntry(self, index, reactant1, product1, kinetics, reactant2=None, reactant3=None, product2=None, product3=None, degeneracy=1, label='', duplicate=False, reference=None, referenceType='', shortDesc='', longDesc='', history=None):
        
        reactants = [Species(label=reactant1.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant1)])]
        if reactant2 is not None: reactants.append(Species(label=reactant2.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant2)]))
        if reactant3 is not None: reactants.append(Species(label=reactant3.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(reactant3)]))

        products = [Species(label=product1.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product1)])]
        if product2 is not None: products.append(Species(label=product2.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product2)]))
        if product3 is not None: products.append(Species(label=product3.strip().splitlines()[0].strip(), molecule=[Molecule().fromAdjacencyList(product3)]))
        
        comment = "Reaction from {0}".format(self.label)
        if shortDesc.strip(): 
            comment += "{0!s}\n".format(shortdesc.strip())
        if longDesc.strip():
            comment += str(re.sub('\s*\n\s*','\n',longDesc))
        kinetics.comment = comment.strip()
        
        self.entries[index] = Entry(
            index = index,
            label = label,
            item = Reaction(reactants=reactants, products=products, degeneracy=degeneracy, duplicate=duplicate),
            data = kinetics,
            reference = reference,
            referenceType = referenceType,
            shortDesc = shortDesc,
            longDesc = longDesc.strip(),
            history = history or [],
        )

    def saveEntry(self, f, entry):
        """
        Write the given `entry` in the kinetics library to the file object `f`.
        """
        return saveEntry(f, entry)

    def loadOld(self, path):
        """
        Load an old-style RMG kinetics library from the location `path`.
        """
        path = os.path.abspath(path)

        self.loadOldDictionary(os.path.join(path,'species.txt'), pattern=False)
        species = dict([(label, Species(label=label, molecule=[entry.item])) for label, entry in self.entries.iteritems()])
        
        reactions = []
        reactions.extend(self.__loadOldReactions(os.path.join(path,'reactions.txt'), species))
        if os.path.exists(os.path.join(path,'pdepreactions.txt')):
            reactions.extend(self.__loadOldReactions(os.path.join(path,'pdepreactions.txt'), species))

        self.entries = {}
        for index, reaction in enumerate(reactions):
            entry = Entry(
                index = index+1,
                item = reaction,
                data = reaction.kinetics,
            )
            entry.longDesc = reaction.kinetics.comment
            reaction.kinetics.comment = ''
            self.entries[index+1] = entry
            reaction.kinetics = None
        
        self.checkForDuplicates()
        self.convertDuplicatesToMulti()

    def __loadOldReactions(self, path, species):
        """
        Load an old-style reaction library from `path`. This algorithm can
        handle both the pressure-independent and pressure-dependent reaction
        files. If the pressure-dependent file is read, the extra pressure-
        dependent kinetics information is ignored unless the kinetics database
        is a seed mechanism.
        """
        reactions = []
        
        # Process the reactions or pdepreactions file
        try:
            inUnitSection = False; inReactionSection = False
            Aunits = []; Eunits = ''
            reaction = None; kinetics = None
            next_reaction_comment = ''

            fdict = open(path, 'r')
            for line in fdict:
                line, comment = splitLineAndComment(line)
                line = line.strip()
                if len(line) == 0:
                    comment = comment.strip()
                    # collect all comment lines and assume they're for the following reaction 
                    next_reaction_comment += comment + '\n'
                    continue
                else: # line is not empty
                    if inUnitSection:
                        if 'A:' in line or 'E:' in line:
                            units = line.split()[1]
                            if 'A:' in line:
                                Aunits0 = units.split('/') # Assume this is a 3-tuple: moles or molecules, volume, time
                                Aunits0[1] = Aunits0[1][0:-1] # Remove '3' from e.g. 'm3' or 'cm3'; this is assumed
                                Aunits = [
                                    '',                                                         # Zeroth-order
                                    '{0}^-1'.format(Aunits0[2]),                                     # First-order
                                    '{0}^3/({1}*{2})'.format(Aunits0[1], Aunits0[0], Aunits0[2]),      # Second-order
                                    '{0}^6/({1}^2*{2})'.format(Aunits0[1], Aunits0[0], Aunits0[2]),    # Third-order
                                ]
                            elif 'E:' in line:
                                Eunits = units
                    elif inReactionSection:
                        if '=' in line:
                            # This line contains a reaction equation, (high-P) Arrhenius parameters, and uncertainties

                            # Strip out "(+M)" from line
                            line = line.replace("(+M)", "")

                            items = line.split()

                            # Find the reaction arrow
                            for arrow in ['<=>', '=>', '=', '->']:
                                if arrow in items:
                                    arrowIndex = items.index(arrow)
                                    break

                            # Find the start of the data
                            try:
                                temp = float(items[-6])
                                dataIndex = -6
                            except ValueError:
                                dataIndex = -3

                            # Find the reactant and product items
                            reactantItems = []
                            for item in items[0:arrowIndex]:
                                if item != '+':
                                    for i in item.split('+'):
                                        if i != '' and i != 'M': reactantItems.append(i)
                            productItems = []
                            for item in items[arrowIndex+1:dataIndex]:
                                if item != '+':
                                    for i in item.split('+'):
                                        if i != '' and i != 'M': productItems.append(i)
                            
                            reactants = []; products = []
                            for item in reactantItems:
                                try:
                                    reactants.append(species[item])
                                except KeyError:
                                    raise DatabaseError('Reactant {0} not found in species dictionary.'.format(item))
                            for item in productItems:
                                try:
                                    products.append(species[item])
                                except KeyError:
                                    raise DatabaseError('Product {0} not found in species dictionary.'.format(item))

                            if dataIndex == -6:
                                A, n, Ea, dA, dn, dEa = items[-6:]
                                A = float(A)
                            else:
                                A, n, Ea = items[-3:]
                                dA = '0'; dn = '0'; dEa = '0'
                            
                            A = float(A)
                            if dA[0] == '*':
                                A = Quantity(A,Aunits[len(reactants)],'*|/',float(dA[1:]))
                            else:
                                dA = float(dA)
                                if dA != 0:
                                    A = Quantity(A,Aunits[len(reactants)],'+|-',dA)
                                else:
                                    A = Quantity(A,Aunits[len(reactants)])

                            n = float(n); dn = float(dn)
                            if dn != 0:
                                n = Quantity(n,'','+|-',dn)
                            else:
                                n = Quantity(n,'')

                            Ea = float(Ea); dEa = float(dEa)
                            if dEa != 0:
                                Ea = Quantity(Ea,Eunits,'+|-',dEa)
                            else:
                                Ea = Quantity(Ea,Eunits)

                            kinetics = Arrhenius(A=A, n=n, Ea=Ea, T0=(1.0,"K"))

                            reaction = Reaction(
                                reactants=reactants,
                                products=products,
                                kinetics=kinetics,
                                reversible=(arrow in ['<=>', '=']),
                            )
                            reaction.kinetics.comment = next_reaction_comment
                            next_reaction_comment = ""
                            reactions.append(reaction)

                        elif 'PLOG' in line:
                            # This line contains pressure-dependent Arrhenius parameters in Chemkin format
                            items = line.split('/')
                            P, A, n, Ea = items[1].split()
                            P = float(P)
                            A = Quantity(float(A), Aunits[len(reactants)])
                            n = Quantity(float(n), '')
                            Ea = Quantity(float(Ea), Eunits)
                            arrhenius = Arrhenius(A=A, n=n, Ea=Ea, T0=(1.0,"K"))
                            if not isinstance(kinetics, PDepArrhenius):
                                old_kinetics = kinetics
                                comment = old_kinetics.comment
                                old_kinetics.comment = ''
                                assert isinstance(old_kinetics, Arrhenius)
                                kinetics = PDepArrhenius(pressures=([P],"atm"), arrhenius=[arrhenius], highPlimit=old_kinetics, comment=comment)
                            else:
                                pressures = list(kinetics.pressures.values)
                                pressures.append(P*101325.)
                                kinetics.pressures.values = numpy.array(pressures, numpy.float)
                                kinetics.arrhenius.append(arrhenius)
                            reaction.kinetics = kinetics
                            
                        elif 'LOW' in line:
                            # This line contains low-pressure-limit Arrhenius parameters in Chemkin format

                            # Upgrade the kinetics to a Lindemann if not already done
                            if isinstance(kinetics, ThirdBody):
                                kinetics = Lindemann(arrheniusHigh=kinetics.arrheniusHigh,
                                                     efficiencies=kinetics.efficiencies,
                                                     comment=kinetics.comment)
                                reaction.kinetics = kinetics
                            elif isinstance(kinetics, Arrhenius):
                                kinetics = Lindemann(arrheniusHigh=kinetics, comment=kinetics.comment)
                                kinetics.arrheniusHigh.comment = ''
                                reaction.kinetics = kinetics

                            items = line.split('/')
                            A, n, Ea = items[1].split()
                            A = Quantity(float(A), Aunits[len(reactants)])
                            n = Quantity(float(n), '')
                            Ea = Quantity(float(Ea), Eunits)
                            kinetics.arrheniusLow = Arrhenius(A=A, n=n, Ea=Ea, T0=1.0)

                        elif 'TROE' in line:
                            # This line contains Troe falloff parameters in Chemkin format

                            # Upgrade the kinetics to a Troe if not already done
                            if isinstance(kinetics, Lindemann):
                                kinetics = Troe(arrheniusLow=kinetics.arrheniusLow,
                                                arrheniusHigh=kinetics.arrheniusHigh,
                                                efficiencies=kinetics.efficiencies,
                                                comment=kinetics.comment)
                                reaction.kinetics = kinetics
                            elif isinstance(kinetics, ThirdBody):
                                kinetics = Troe(arrheniusHigh=kinetics.arrheniusHigh,
                                                efficiencies=kinetics.efficiencies,
                                                comment=kinetics.comment)
                                reaction.kinetics = kinetics
                            elif isinstance(kinetics, Arrhenius):
                                kinetics = Troe(arrheniusHigh=kinetics, comment=kinetics.comment)
                                kinetics.arrheniusHigh.comment = ''
                                reaction.kinetics = kinetics

                            items = line.split('/')
                            items = items[1].split()
                            if len(items) == 3:
                                alpha, T3, T1 = items; T2 = None
                            else:
                                alpha, T3, T1, T2 = items

                            kinetics.alpha = Quantity(float(alpha))
                            kinetics.T1 = Quantity(float(T1),"K")
                            if T2 is not None:
                                kinetics.T2 = Quantity(float(T2),"K")
                            else:
                                kinetics.T2 = None
                            kinetics.T3 = Quantity(float(T3),"K")

                        elif 'DUPLICATE' in line or 'DUP' in line:
                            reaction.duplicate = True
                        
                        else:
                            # This line contains collider efficiencies

                            # Upgrade the kinetics to a ThirdBody if not already done
                            if isinstance(kinetics, Arrhenius):
                                kinetics = ThirdBody(arrheniusHigh=kinetics, comment=kinetics.comment)
                                kinetics.arrheniusHigh.comment = ''
                                reaction.kinetics = kinetics

                            items = line.split('/')
                            for spec, eff in zip(items[0::2], items[1::2]):
                                spec = str(spec).strip()

                                # In old database, N2, He, Ne, and Ar were treated as special "bath gas" species
                                # These bath gas species were not required to be in the species dictionary
                                # The new database removes this special case, and requires all colliders to be explicitly defined
                                # This is hardcoding to handle these special colliders
                                if spec.upper() in ['N2', 'HE', 'AR', 'NE'] and spec not in species:
                                    if spec.upper() == 'N2':
                                        species[spec] = Species(label='N2', molecule=[Molecule().fromSMILES('N#N')])
                                    elif spec.upper() == 'HE':
                                        species[spec] = Species(label='He', molecule=[Molecule().fromAdjacencyList('1 He 0')])
                                    elif spec.upper() == 'AR':
                                        species[spec] = Species(label='Ar', molecule=[Molecule().fromAdjacencyList('1 Ar 0')])
                                    elif spec.upper() == 'NE':
                                        species[spec] = Species(label='Ne', molecule=[Molecule().fromAdjacencyList('1 Ne 0')])
                                
                                if spec not in species:
                                    logging.warning('Collider {0} for reaction {1} not found in species dictionary.'.format(spec, reaction))
                                else:
                                    kinetics.efficiencies[species[spec].molecule[0]] = float(eff)

                    if 'Unit:' in line:
                        inUnitSection = True; inReactionSection = False
                    elif 'Reactions:' in line:
                        inUnitSection = False; inReactionSection = True
                        
        except (DatabaseError, InvalidAdjacencyListError), e:
            logging.exception(str(e))
            raise
        except IOError, e:
            logging.exception('Database dictionary file "' + e.filename + '" not found.')
            raise
        finally:
            if fdict: fdict.close()

        return reactions
    
    def saveOld(self, path):
        """
        Save an old-style reaction library to `path`. This creates files named
        ``species.txt``, ``reactions.txt``, and ``pdepreactions.txt`` in the
        given directory; these contain the species dictionary, high-pressure
        limit reactions and kinetics, and pressure-dependent reactions and
        kinetics, respectively.
        """
        try:
            os.makedirs(path)
        except OSError:
            pass
        
        def writeArrhenius(f, arrhenius):
            f.write('{0:10.3e} {1:9.3f} {2:10.2f}     {3}{4:g} {5:g} {6:g}\n'.format(
                arrhenius.A.value,
                arrhenius.n.value,
                arrhenius.Ea.value / 4.184,
                '*' if arrhenius.A.isUncertaintyMultiplicative() else '',
                arrhenius.A.uncertainty,
                arrhenius.n.uncertainty,
                arrhenius.Ea.uncertainty / 4.184,
            ))
        
        # Gather all of the species used in this kinetics library
        speciesDict = self.getSpecies()
        # Also include colliders in the above
        for entry in self.entries.values():
            if isinstance(entry.data, ThirdBody):
                for molecule in entry.data.efficiencies:
                    formula = molecule.getFormula()
                    if formula in ['He', 'Ar', 'N2', 'Ne']:
                        pass
                    else:
                        found = False
                        for species in speciesDict.values():
                            for mol in species.molecule:
                                if mol.isIsomorphic(molecule):
                                    found = True
                                    break
                        if not found:
                            speciesDict[formula] = Species(label=formula, molecule=[molecule])
        
        entries = self.entries.values()
        entries.sort(key=lambda x: x.index)
        
        # Save the species dictionary
        speciesList = speciesDict.values()
        speciesList.sort(key=lambda x: x.label)
        f = open(os.path.join(path, 'species.txt'), 'w')
        for species in speciesList:
            f.write(species.molecule[0].toAdjacencyList(label=species.label, removeH=True) + "\n")
        f.close()
        
        # Save the high-pressure limit reactions
        # Currently only Arrhenius kinetics are allowed
        f = open(os.path.join(path, 'reactions.txt'), 'w')
        f.write('Unit:\n')
        f.write('A: mol/m3/s\n')
        f.write('E: cal/mol\n\n')
        f.write('Reactions:\n')
        for entry in entries:
            kinetics = entry.data
            if not kinetics.isPressureDependent():
                # Write reaction equation
                f.write('{0:<48}'.format(entry.item))
                # Write kinetics
                if isinstance(kinetics, Arrhenius):
                    writeArrhenius(f, kinetics)
                else:
                    raise DatabaseError('Unexpected kinetics type "{0}" encountered while saving old kinetics library (reactions.txt).'.format(kinetics.__class__))
                # Mark as duplicate if needed
                if entry.item.duplicate:
                    f.write(' DUPLICATE\n')
        f.close()
        
        # Save the pressure-dependent reactions
        # Currently only ThirdBody, Lindemann, Troe, and PDepArrhenius kinetics are allowed
        f = open(os.path.join(path, 'pdepreactions.txt'), 'w')
        f.write('Unit:\n')
        f.write('A: mol/m3/s\n')
        f.write('E: cal/mol\n\n')
        f.write('Reactions:\n')
        for entry in entries:
            kinetics = entry.data
            if entry.data.isPressureDependent():
                # Write reaction equation
                equation = str(entry.item)
                index = equation.find('<=>')
                if isinstance(kinetics, ThirdBody) and not isinstance(kinetics, Lindemann):
                    equation = '{0}+ M {1} + M'.format(equation[0:index], equation[index:])
                elif isinstance(kinetics, PDepArrhenius):
                    pass
                else:
                    equation = '{0}(+M) {1} (+M)'.format(equation[0:index], equation[index:]) 
                f.write('{0:<48}'.format(equation))
                # Write kinetics
                if isinstance(kinetics, ThirdBody):
                    writeArrhenius(f, kinetics.arrheniusHigh)
                    if len(kinetics.efficiencies) > 0:
                        for molecule, efficiency in kinetics.efficiencies.iteritems():
                            for spec in speciesDict.values():
                                if molecule in spec.molecule:
                                    f.write('{0}/{1:.2f}/ '.format(spec.label, efficiency)) 
                            else:
                                f.write('{0}/{1:.2f}/ '.format(molecule.getFormula().upper(), efficiency)) 
                        f.write('\n')
                    if isinstance(kinetics, Lindemann):
                        f.write('     LOW  /  {0:10.3e} {1:9.3f} {2:10.2f}/\n'.format(
                            kinetics.arrheniusLow.A.value,
                            kinetics.arrheniusLow.n.value,
                            kinetics.arrheniusLow.Ea.value / 4.184,
                        ))
                    if isinstance(kinetics, Troe):
                        if kinetics.T2 is not None:
                            f.write('     TROE /  {0:10.4f} {1:10.2g} {2:10.2g} {3:10.2g}/\n'.format(
                                kinetics.alpha.value,
                                kinetics.T3.value,
                                kinetics.T1.value,
                                kinetics.T2.value,
                            ))
                        else:
                            f.write('     TROE /  {0:10.4f} {1:10.2g} {2:10.2g}/\n'.format(
                                kinetics.alpha.value,
                                kinetics.T3.value,
                                kinetics.T1.value,
                            ))
                        
                elif isinstance(kinetics, PDepArrhenius):
                    for pressure, arrhenius in zip(kinetics.pressures.values, kinetics.arrhenius):
                        f.write('     PLOG /  {0:10g} {1:10.3e} {2:9.3f} {3:10.2f} /\n'.format(
                            pressure / 1e5,
                            arrhenius.A.value,
                            arrhenius.n.value,
                            arrhenius.Ea.value / 4.184,
                        ))
                else:
                    raise DatabaseError('Unexpected kinetics type "{0}" encountered while saving old kinetics library (reactions.txt).'.format(kinetics.__class__))
                # Mark as duplicate if needed
                if entry.item.duplicate:
                    f.write(' DUPLICATE\n')
                f.write('\n')
        f.close()
    
################################################################################

class KineticsGroups(Database):
    """
    A class for working with an RMG kinetics family group additivity values. 
    """

    def __init__(self, entries=None, top=None, label='', name='', shortDesc='', longDesc='', forwardTemplate=None, forwardRecipe=None, reverseTemplate=None, reverseRecipe=None, forbidden=None):
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
                else:
                    logging.warning("Couldn't find match for {0} in {1}".format(entry,atomList))
                    logging.warning(reactant.toAdjacencyList())

        # Get fresh templates (with duplicate nodes back in)
        forwardTemplate = self.top[:]
        if self.label.lower().startswith('r_recombination'):
            forwardTemplate.append(forwardTemplate[0])

        # Check that we were able to match the template.
        # template is a list of the actual matched nodes
        # forwardTemplate is a list of the top level nodes that should be matched
        if len(template) != len(forwardTemplate):
            logging.warning('Unable to find matching template for reaction {0} in reaction family {1}'.format(str(reaction), str(self)) )
            logging.warning(" Trying to match " + str(forwardTemplate))
            logging.warning(" Matched "+str(template))
            print str(self), template, forwardTemplate
            for reactant in reaction.reactants:
                print reactant.toAdjacencyList() + '\n'
            for product in reaction.products:
                print product.toAdjacencyList() + '\n'
            raise UndeterminableKineticsError(reaction)

        return template

    def getKineticsForTemplate(self, template, referenceKinetics, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template`.
        """

        # Start with the generic kinetics of the top-level nodes
        kinetics = referenceKinetics
        
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
            kinetics.kdata.values *= degeneracy
        elif isinstance(kinetics, Arrhenius):
            kinetics.A.value *= degeneracy
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
            if len(kinetics1.Tdata.values) != len(kinetics2.Tdata.values) or any([T1 != T2 for T1, T2 in zip(kinetics1.Tdata.values, kinetics2.Tdata.values)]):
                raise KineticsError('Cannot add these KineticsData objects due to their having different temperature points.')
            kinetics = KineticsData(
                Tdata = (kinetics1.Tdata.values, kinetics2.Tdata.units),
                kdata = (kinetics1.kdata.values * kinetics2.kdata.values, kinetics1.kdata.units),
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
            kinetics.Tmin = kinetics1.Tmin if kinetics1.Tmin.value > kinetics2.Tmin.value else kinetics2.Tmin
        elif kinetics1.Tmin is not None and kinetics2.Tmin is None:
            kinetics.Tmin = kinetics1.Tmin
        elif kinetics1.Tmin is None and kinetics2.Tmin is not None:
            kinetics.Tmin = kinetics2.Tmin
        
        if kinetics1.Tmax is not None and kinetics2.Tmax is not None:
            kinetics.Tmax = kinetics1.Tmax if kinetics1.Tmax.value < kinetics2.Tmax.value else kinetics2.Tmax
        elif kinetics1.Tmax is not None and kinetics2.Tmax is None:
            kinetics.Tmax = kinetics1.Tmax
        elif kinetics1.Tmax is None and kinetics2.Tmax is not None:
            kinetics.Tmax = kinetics2.Tmax
        
        if kinetics1.Pmin is not None and kinetics2.Pmin is not None:
            kinetics.Pmin = kinetics1.Pmin if kinetics1.Pmin.value > kinetics2.Pmin.value else kinetics2.Pmin
        elif kinetics1.Pmin is not None and kinetics2.Pmin is None:
            kinetics.Pmin = kinetics1.Pmin
        elif kinetics1.Pmin is None and kinetics2.Pmin is not None:
            kinetics.Pmin = kinetics2.Pmin
        
        if kinetics1.Pmax is not None and kinetics2.Pmax is not None:
            kinetics.Pmax = kinetics1.Pmax if kinetics1.Pmax.value < kinetics2.Pmax.value else kinetics2.Pmax
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

    def __init__(self, entries=None, top=None, label='', name='', shortDesc='', longDesc='', forwardTemplate=None, forwardRecipe=None, reverseTemplate=None, reverseRecipe=None, forbidden=None):
        Database.__init__(self, entries, top, label, name, shortDesc, longDesc)
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

        self.groups.loadOldDictionary(os.path.join(path, 'dictionary.txt'), pattern=True)
        self.groups.loadOldTree(os.path.join(path, 'tree.txt'))
        # The old kinetics groups use rate rules (not group additivity values),
        # so we can't load the old rateLibrary.txt

        # Load the reaction recipe
        self.loadOldTemplate(os.path.join(path, 'reactionAdjList.txt'))
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
        if os.path.exists(os.path.join(path, 'forbiddenGroups.txt')):
            self.forbidden = ForbiddenStructures().loadOld(os.path.join(path, 'forbiddenGroups.txt'))
            
        entries = self.groups.top[:]
        for entry in self.groups.top:
            entries.extend(self.groups.descendants(entry))
        for index, entry in enumerate(entries):
            entry.index = index + 1
            
        self.rules = KineticsDepository(label='{0}/rules'.format(self.label))
        self.rules.loadOldRateRules(path, self)
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
            ftemp.write('reverse: {0}_reverse\n'.format(self.label))
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
                for f in files:
                    if not f.endswith('.py'): continue
                    name = f.strip('.py')
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
        f.write('shortDesc = "{0}"\n'.format(self.shortDesc))
        f.write('longDesc = """\n')
        f.write(self.longDesc)
        f.write('\n"""\n\n')

        # Write the template
        f.write('template(reactants=[{0}], products=[{1}], ownReverse={2})\n\n'.format(
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.reactants]),
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forwardTemplate.products]),
            self.ownReverse))

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

    def __generateProductStructures(self, reactantStructures, maps, forward):
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

        # If there are two structures and they are the same, then make a copy
        # of the second one and adjust the second map to point to its atoms
        # This is for the case where A + A --> products
        if len(reactantStructures) == 2 and reactantStructures[0] == reactantStructures[1]:
            reactantStructures[1] = reactantStructures[1].copy(deep=True)
            newMap = {}
            for reactantAtom, templateAtom in maps[1].iteritems():
                index = reactantStructures[0].atoms.index(reactantAtom)
                newMap[reactantStructures[1].atoms[index]] = templateAtom
            maps[1] = newMap

        # Tag atoms with labels
        for m in maps:
            for reactantAtom, templateAtom in m.iteritems():
                reactantAtom.label = templateAtom.label

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

        # Check that reactant and product structures are allowed in this family
        # If not, then stop
        if self.forbidden is not None:
            for struct in reactantStructures:
                if self.forbidden.isMoleculeForbidden(struct): raise ForbiddenStructureException()
            for struct in productStructures:
                if self.forbidden.isMoleculeForbidden(struct): raise ForbiddenStructureException()

        # Also check the global forbiddenStructures
        from rmgpy.data.rmg import database
        for struct in reactantStructures:
            if database.forbiddenStructures.isMoleculeForbidden(struct): raise ForbiddenStructureException()
        for struct in productStructures:
            if database.forbiddenStructures.isMoleculeForbidden(struct): raise ForbiddenStructureException()

        return productStructures

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

        # If forbidden structures are defined, make sure the products are not forbidden
        if self.forbidden:
            for product in products:
                if self.forbidden.isMoleculeForbidden(product):
                    return None
        # Also check the global forbiddenStructures
        from rmgpy.data.rmg import database
        for product in products:
            if database.forbiddenStructures.isMoleculeForbidden(product): return None

        # We need to save the reactant and product structures with atom labels so
        # we can generate the kinetics
        # We make copies so the structures aren't trampled on by later actions
        reactants = [reactant.copy(deep=True) for reactant in reactants]
        products = [product.copy(deep=True) for product in products]
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
                ismatch, map = reactant.findSubgraphIsomorphisms(child_structure)
                if ismatch: mappings.extend(map)
            return len(mappings) > 0, mappings
        elif isinstance(struct, Group):
            return reactant.findSubgraphIsomorphisms(struct)

    def generateReactions(self, reactants):
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
        reactions = self.__generateReactions(reactants, forward=True)
        for rxn in reactions:
            reaction = TemplateReaction(
                reactants = rxn.reactants[:],
                products = rxn.products[:],
                degeneracy = rxn.degeneracy,
                thirdBody = rxn.thirdBody,
                reversible = rxn.reversible,
                family = self,
            )
            reactionList.append(reaction)
        
        reverseReactions = []
        if self.ownReverse:
            # for each reaction, make its reverse reaction and store in a 'reverse' attribute
            for rxn in reactionList:
                reactions = self.__generateReactions(rxn.products, forward=True)
                reactions = filterReactions(rxn.products, rxn.reactants, reactions)
                assert len(reactions) == 1, "Expecting one matching reverse reaction, not {0}. Forward reaction {1!s} : {1!r}".format(len(reactions), rxn)
                reaction = reactions[0]
                reaction = TemplateReaction(
                    reactants = reaction.reactants[:],
                    products = reaction.products[:],
                    degeneracy = reaction.degeneracy,
                    thirdBody = reaction.thirdBody,
                    reversible = reaction.reversible,
                    family = self,
                )
                rxn.reverse = reaction
                reverseReactions.append(reaction)
            
        else: # family is not ownReverse
            # Reverse direction (the direction in which kinetics is not defined)
            reactions = self.__generateReactions(reactants, forward=False)
            for rxn in reactions:
                reaction = TemplateReaction(
                    reactants = rxn.products[:],
                    products = rxn.reactants[:],
                    degeneracy = rxn.degeneracy,
                    thirdBody = rxn.thirdBody,
                    reversible = rxn.reversible,
                    family = self,
                )
                reactionList.append(reaction)

        # Also store the reaction template (useful so we can easily get the kinetics later)
        for reaction in reactionList:
            reaction.template = self.getReactionTemplate(reaction)
            if hasattr(reaction,'reverse'):
                reaction.reverse.template = self.getReactionTemplate(reaction.reverse)
        
        # Return the reactions as containing Species objects, not Molecule objects
        for reaction in reactionList:
            reaction.reactants = [Species(label=reactant.toSMILES(), molecule=[reactant]) for reactant in reaction.reactants]
            reaction.products = [Species(label=product.toSMILES(), molecule=[product]) for product in reaction.products]

        return reactionList
    
    def __generateReactions(self, reactants, forward=True):
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

        # Reactants must have explicit hydrogen atoms
        implicitH = []
        for reactant in reactants:
            implicit = []
            for isomer in reactant:
                implicit.append(isomer.implicitHydrogens)
                isomer.makeHydrogensExplicit()
            implicitH.append(implicit)

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

                ismatch, mappings = self.__matchReactantToTemplate(molecule, template.reactants[0])
                if ismatch:
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

                    # Reactants stored as A + B
                    ismatchA, mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[0])
                    ismatchB, mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[1])

                    # Iterate over each pair of matches (A, B)
                    if ismatchA and ismatchB:
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
                        ismatchA, mappingsA = self.__matchReactantToTemplate(moleculeA, template.reactants[1])
                        ismatchB, mappingsB = self.__matchReactantToTemplate(moleculeB, template.reactants[0])

                        # Iterate over each pair of matches (A, B)
                        if ismatchA and ismatchB:
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
        if self.label.lower().startswith('r_recombination'):
            for rxn in rxnList:
                rxn.degeneracy /= 2

        # Restore implicit hydrogen atoms if necessary
        for reactant, implicit in zip(reactants, implicitH):
            for isomer, H in zip(reactant, implicit):
                if H: isomer.makeHydrogensImplicit()

        # This reaction list has only checked for duplicates within itself, not
        # with the global list of reactions
        return rxnList

    def getReactionTemplate(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        return self.groups.getReactionTemplate(reaction)

    def getKineticsForTemplate(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template`.
        """
        # Start with the generic kinetics of the top-level nodes
        kinetics = None
        for entry in self.forwardTemplate.reactants:
            if kinetics is None and entry.data is not None:
                kinetics = deepcopy(entry.data)
        # Now add in more specific corrections if possible
        return self.groups.getKineticsForTemplate(template, kinetics, degeneracy)

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
                    kinetics.A.value *= degeneracy
                    kinetics.comment += "Matched rule {0} {1} in {2}".format(entry.index, entry.label, depository.label)
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
    
    def getKinetics(self, reaction, template, degeneracy=1, returnAllKinetics=True):
        """
        Return the kinetics for the given `reaction` by searching the various
        depositories as well as generating a group additivity estimate. Unlike
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
            for kinetics, entry, isForward in self.getKineticsFromDepository(depository, reaction, template, degeneracy):
                if not returnAllKinetics:
                    return kinetics, depository, entry, isForward
                kineticsList.append([kinetics, depository, entry, isForward])
        # Also generate a group additivity estimate
        kinetics = self.getKineticsForTemplate(template, degeneracy)
        if kinetics:
            if not returnAllKinetics:
                return kinetics, None, None, True
            kineticsList.append([kinetics, None, None, True])
        
        if not returnAllKinetics:
            raise UndeterminableKineticsError(reaction)
        
        return kineticsList

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
        
###########
class KineticsDatabase:
    """
    A class for working with the RMG kinetics database.
    """

    def __init__(self):
        self.families = {}
        self.libraries = {}
        self.libraryOrder = []
        self.local_context = {
            'KineticsData': KineticsData,
            'Arrhenius': Arrhenius,
            'ArrheniusEP': ArrheniusEP,
            'MultiKinetics': MultiKinetics,
            'PDepArrhenius': PDepArrhenius,
            'Chebyshev': Chebyshev,
            'ThirdBody': ThirdBody,
            'Lindemann': Lindemann,
            'Troe': Troe,
        }
        self.global_context = {}

    def load(self, path, libraries=None, depositories=None):
        """
        Load the kinetics database from the given `path` on disk, where `path`
        points to the top-level folder of the families database.
        """
        self.loadFamilies(os.path.join(path, 'families'), depositories)
        self.loadLibraries(os.path.join(path, 'libraries'), libraries)
        
    def loadFamilies(self, path, depositories=None):
        """
        Load the kinetics families from the given `path` on disk, where `path`
        points to the top-level folder of the kinetics families.
        """
        self.families = {}
        logging.info('Loading kinetics families from {0}'.format(path))
        for (root, dirs, files) in os.walk(os.path.join(path)):
            for d in dirs:
                familyPath = os.path.join(root, d)
                family = KineticsFamily(label=d)
                family.load(familyPath, self.local_context, self.global_context, depositoryLabels=depositories)
                self.families[d] = family

    def loadLibraries(self, path, libraries=None):
        """
        Load the listed kinetics libraries from the given `path` on disk.
        
        Loads them all if `libraries` list is not specified or `None`.
        The `path` points to the folder of kinetics libraries in the database,
        and the libraries should be in files like :file:`<path>/<library>.py`.
        """
        self.libraries = {}; self.libraryOrder = []
        
        if libraries is not None:
            for library_name in libraries:
                library_file = os.path.join(path, library_name+'.py')
                if os.path.exists(library_file):
                    logging.info('Loading kinetics library {0} from {1}...'.format(library_name, library_file))
                    library = KineticsLibrary(label=library_name)
                    library.load(library_file, self.local_context, self.global_context)
                    self.libraries[library.label] = library
                    self.libraryOrder.append(library.label)
                else:
                    raise IOError("Couldn't find kinetics library {0}".format(library_file))
            assert (self.libraryOrder == libraries)
        else:# load all the libraries you can find
            for (root, dirs, files) in os.walk(os.path.join(path)):
                for f in files:
                    name, ext = os.path.splitext(f)
                    if ext.lower() == '.py':
                        library_file = os.path.join(root, f)
                        label=library_file[len(path)+1:-3]
                        logging.info('Loading kinetics library {0} from {1}...'.format(label, library_file))
                        library = KineticsLibrary(label=label)
                        library.load(library_file, self.local_context, self.global_context)
                        self.libraries[library.label] = library
                        self.libraryOrder.append(library.label)

    def save(self, path):
        """
        Save the kinetics database to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics database.
        """
        path = os.path.abspath(path)
        if not os.path.exists(path): os.mkdir(path)
        self.saveFamilies(os.path.join(path, 'families'))
        self.saveLibraries(os.path.join(path, 'libraries'))

    def saveFamilies(self, path):
        """
        Save the kinetics families to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics families.
        """
        if not os.path.exists(path): os.mkdir(path)
        for label, family in self.families.iteritems():
            familyPath = os.path.join(path, label)
            if not os.path.exists(familyPath): os.mkdir(familyPath)
            family.save(familyPath)

    def saveLibraries(self, path):
        """
        Save the kinetics libraries to the given `path` on disk, where `path`
        points to the top-level folder of the kinetics libraries.
        """
        for label, library in self.libraries.iteritems():
            folders = label.split(os.sep)
            try:
                os.makedirs(os.path.join(path, *folders[:-1]))
            except OSError:
                pass
            library.save(os.path.join(path, '{0}.py'.format(label)))

    def loadOld(self, path):
        """
        Load the old RMG kinetics database from the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        self.families = {}
        self.libraries = {}
        
        librariesPath = os.path.join(path, 'kinetics_libraries')
        for (root, dirs, files) in os.walk(os.path.join(path, 'kinetics_libraries')):
            if os.path.exists(os.path.join(root, 'species.txt')) and os.path.exists(os.path.join(root, 'reactions.txt')):
                library = KineticsLibrary(label=root[len(librariesPath)+1:], name=root[len(librariesPath)+1:])
                library.loadOld(root)
                self.libraries[library.label] = library
                
        for (root, dirs, files) in os.walk(os.path.join(path, 'kinetics_groups')):
            if os.path.exists(os.path.join(root, 'dictionary.txt')) and os.path.exists(os.path.join(root, 'rateLibrary.txt')):
                label = os.path.split(root)[1]
                family = KineticsFamily(label=label)
                family.loadOld(root)
                self.families[family.label] = family

        return self

    def saveOld(self, path):
        """
        Save the old RMG kinetics database to the given `path` on disk, where
        `path` points to the top-level folder of the old RMG database.
        """
        librariesPath = os.path.join(path, 'kinetics_libraries')
        if not os.path.exists(librariesPath): os.mkdir(librariesPath)
        for library in self.libraries.values():
            libraryPath = os.path.join(librariesPath, library.label)
            library.saveOld(libraryPath)

        groupsPath = os.path.join(path, 'kinetics_groups')
        if not os.path.exists(groupsPath): os.mkdir(groupsPath)
        for label, family in self.families.iteritems():
            groupPath = os.path.join(groupsPath, label)
            family.saveOld(groupPath)
    
    def generateReactions(self, reactants, products=None):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository, libraries, and groups, in that order.
        """
        reactionList = []
        reactionList.extend(self.generateReactionsFromLibraries(reactants, products))
        reactionList.extend(self.generateReactionsFromFamilies(reactants, products))
        return reactionList

    def generateReactionsFromLibraries(self, reactants, products):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository.
        """
        reactionList = []
        for label in self.libraryOrder:
            reactionList.extend(self.generateReactionsFromLibrary(reactants, products, self.libraries[label]))
        return reactionList

    def generateReactionsFromLibrary(self, reactants, products, library):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        searches the depository.
        """
        reactionList = []
        for entry in library.entries.values():
            if entry.item.matchesMolecules(reactants):
                reaction = LibraryReaction(
                    reactants = entry.item.reactants[:],
                    products = entry.item.products[:],
                    degeneracy = entry.item.degeneracy,
                    thirdBody = entry.item.thirdBody,
                    reversible = entry.item.reversible,
                    kinetics = deepcopy(entry.data),
                    library = library,
                    entry = entry,
                )
                reactionList.append(reaction)
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        return reactionList

    def generateReactionsFromFamilies(self, reactants, products, only_families=None):
        """
        Generate all reactions between the provided list of one or two
        `reactants`, which should be :class:`Molecule` objects. This method
        applies the reaction family.
        If `only_families` is a list of strings, only families with those labels
        are used.
        """
        reactionList = []
        for label, family in self.families.iteritems():
            if only_families is None or label in only_families:
                reactionList.extend(family.generateReactions(reactants))
        if products:
            reactionList = filterReactions(reactants, products, reactionList)
        return reactionList

    def getForwardReactionForFamilyEntry(self, entry, family, thermoDatabase):
        """
        For a given `entry` for a reaction of the given reaction `family` (the
        string label of the family), return the reaction with kinetics and
        degeneracy for the "forward" direction as defined by the reaction 
        family. For families that are their own reverse, the direction the
        kinetics is given in will be preserved. If the entry contains 
        functional groups for the reactants, assume that it is given in the 
        forward direction and do nothing. Returns the reaction in the direction
        consistent with the reaction family template, and the matching template.
        Note that the returned reaction will have its kinetics and degeneracy
        set appropriately.
        
        In order to reverse the reactions that are given in the reverse of the
        direction the family is defined, we need to compute the thermodynamics
        of the reactants and products. For this reason you must also pass
        the `thermoDatabase` to use to generate the thermo data.
        """

        def generateThermoData(species, thermoDatabase):
            thermoData = [thermoDatabase.getThermoData(molecule) for molecule in species.molecule]
            thermoData.sort(key=lambda x: x.getEnthalpy(298))
            return thermoData[0]
        
        def matchSpeciesToMolecules(species, molecules):
            if len(species) == len(molecules) == 1:
                return species[0].isIsomorphic(molecules[0])
            elif len(species) == len(molecules) == 2:
                if species[0].isIsomorphic(molecules[0]) and species[1].isIsomorphic(molecules[1]):
                    return True
                elif species[0].isIsomorphic(molecules[1]) and species[1].isIsomorphic(molecules[0]):
                    return True
            return False

        reaction = None; template = None

        # Get the indicated reaction family
        try:
            groups = self.families[family].groups
        except KeyError:
            raise ValueError('Invalid value "{0}" for family parameter.'.format(family))

        if all([(isinstance(reactant, Group) or isinstance(reactant, LogicNode)) for reactant in entry.item.reactants]):
            # The entry is a rate rule, containing functional groups only
            # By convention, these are always given in the forward direction and
            # have kinetics defined on a per-site basis
            reaction = Reaction(
                reactants = entry.item.reactants[:],
                products = [],
                kinetics = entry.data,
                degeneracy = 1,
            )
            template = [groups.entries[label] for label in entry.label.split(';')]

        elif (all([isinstance(reactant, Molecule) for reactant in entry.item.reactants]) and
            all([isinstance(product, Molecule) for product in entry.item.products])):
            # The entry is a real reaction, containing molecules
            # These could be defined for either the forward or reverse direction
            # and could have a reaction-path degeneracy

            reaction = Reaction(reactants=[], products=[])
            for molecule in entry.item.reactants:
                molecule.makeHydrogensExplicit()
                reactant = Species(molecule=[molecule], label=molecule.toSMILES())
                reactant.generateResonanceIsomers()
                reactant.thermo = generateThermoData(reactant, thermoDatabase)
                reaction.reactants.append(reactant)
            for molecule in entry.item.products:
                molecule.makeHydrogensExplicit()
                product = Species(molecule=[molecule], label=molecule.toSMILES())
                product.generateResonanceIsomers()
                product.thermo = generateThermoData(product, thermoDatabase)
                reaction.products.append(product)

            # Generate all possible reactions involving the reactant species
            generatedReactions = self.generateReactionsFromFamilies([reactant.molecule for reactant in reaction.reactants], only_families=[family])

            # Remove from that set any reactions that don't produce the desired reactants and products
            forward = []; reverse = []
            for rxn in generatedReactions:
                if matchSpeciesToMolecules(reaction.reactants, rxn.reactants) and matchSpeciesToMolecules(reaction.products, rxn.products):
                    forward.append(rxn)
                if matchSpeciesToMolecules(reaction.reactants, rxn.products) and matchSpeciesToMolecules(reaction.products, rxn.reactants):
                    reverse.append(rxn)

            # We should now know whether the reaction is given in the forward or
            # reverse direction
            if len(forward) == 1 and len(reverse) == 0:
                # The reaction is in the forward direction, so use as-is
                reaction = forward[0]
                template = groups.getReactionTemplate(reaction)
                # Don't forget to overwrite the estimated kinetics from the database with the kinetics for this entry
                reaction.kinetics = entry.data
            elif len(reverse) == 1 and len(forward) == 0:
                # The reaction is in the reverse direction
                # First fit Arrhenius kinetics in that direction
                Tdata = 1.0/numpy.arange(0.0005,0.0035,0.0001,numpy.float64)
                kdata = []
                for T in Tdata:
                    kdata.append(entry.data.getRateCoefficient(T) / reaction.getEquilibriumConstant(T))
                kdata = numpy.array(kdata, numpy.float64)
                kunits = 'm^3/(mol*s)' if len(reverse[0].reactants) == 2 else 's^-1'
                kinetics = Arrhenius().fitToData(Tdata, kdata, kunits, T0=1.0)
                # Now flip the direction
                reaction = reverse[0]
                reaction.kinetics = kinetics
                template = groups.getReactionTemplate(reaction)
            elif len(reverse) > 0 and len(forward) > 0:
                print 'FAIL: Multiple reactions found for "%s".' % (entry.label)
            elif len(reverse) == 0 and len(forward) == 0:
                print 'FAIL: No reactions found for "%s".' % (entry.label)
            else:
                print 'FAIL: Unable to estimate kinetics for "%s".' % (entry.label)

        assert reaction is not None
        assert template is not None
        return reaction, template
