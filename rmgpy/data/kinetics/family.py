#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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
import codecs
import itertools
import logging
import multiprocessing as mp
import os.path
import random
import math
import re
import warnings
from collections import OrderedDict
from copy import deepcopy

import numpy as np
from sklearn.model_selection import KFold

from rmgpy import settings
from rmgpy.constraints import fails_species_constraints
from rmgpy.data.base import Database, Entry, LogicNode, LogicOr, ForbiddenStructures, get_all_combinations
from rmgpy.data.kinetics.common import save_entry, find_degenerate_reactions, generate_molecule_combos, \
                                       ensure_independent_atom_ids, check_for_same_reactants
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.data.kinetics.groups import KineticsGroups
from rmgpy.data.kinetics.rules import KineticsRules
from rmgpy.exceptions import ActionError, DatabaseError, InvalidActionError, KekulizationError, KineticsError, \
                             ForbiddenStructureException, UndeterminableKineticsError, AtomTypeError
from rmgpy.kinetics import Arrhenius, SurfaceArrhenius, SurfaceArrheniusBEP, StickingCoefficient, \
                           StickingCoefficientBEP, ArrheniusBM, SurfaceChargeTransfer, ArrheniusChargeTransfer, \
                           ArrheniusChargeTransferBM, KineticsModel, Marcus
from rmgpy.kinetics.uncertainties import RateUncertainty, rank_accuracy_map
from rmgpy.molecule import Bond, GroupBond, Group, Molecule
from rmgpy.molecule.atomtype import ATOMTYPES
from rmgpy.reaction import Reaction, same_species_lists
from rmgpy.species import Species
from rmgpy.tools.uncertainty import KineticParameterUncertainty
from rmgpy.molecule.fragment import Fragment
import rmgpy.constants as constants
from rmgpy.data.solvation import SoluteData, add_solute_data, SoluteTSData, to_soluteTSdata

################################################################################


class TemplateReaction(Reaction):
    """
    A Reaction object generated from a reaction family template. In addition
    to attributes inherited from :class:`Reaction`, this class includes the
    following attributes:

    =============== ========================= =====================================
    Attribute       Type                      Description
    =============== ========================= =====================================
    `family`        ``str``                   The kinetics family that the reaction was created from.
    `estimator`     ``str``                   The name of the kinetic estimator; currently only rate rules is supported.
    `reverse`       :class:`TemplateReaction` The reverse reaction, for families that are their own reverse.
    `is_forward`    ``bool``                  Whether the reaction was generated in the forward direction of the family.
    `labeled_atoms` ``dict``                  Keys are 'reactants' or 'products', values are dictionaries.
                                              Keys in the second level dictionary are template labels (e.g., ``*1``),
                                              values are the respective Atom object instance in the reactants.
    =============== ========================= =====================================
    """

    def __init__(self,
                 index=-1,
                 reactants=None,
                 products=None,
                 specific_collider=None,
                 kinetics=None,
                 reversible=True,
                 transition_state=None,
                 duplicate=False,
                 degeneracy=1,
                 pairs=None,
                 family=None,
                 template=None,
                 estimator=None,
                 reverse=None,
                 is_forward=None,
                 electrons=0,
                 ):
        Reaction.__init__(self,
                          index=index,
                          reactants=reactants,
                          products=products,
                          specific_collider=specific_collider,
                          kinetics=kinetics,
                          reversible=reversible,
                          transition_state=transition_state,
                          duplicate=duplicate,
                          degeneracy=degeneracy,
                          pairs=pairs,
                          is_forward=is_forward,
                          electrons=electrons
                          )
        self.family = family
        self.template = template
        self.estimator = estimator
        self.reverse = reverse
        self.labeled_atoms = {'reactants': dict(), 'products': dict()}

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (TemplateReaction, (self.index,
                                   self.reactants,
                                   self.products,
                                   self.specific_collider,
                                   self.kinetics,
                                   self.reversible,
                                   self.transition_state,
                                   self.duplicate,
                                   self.degeneracy,
                                   self.pairs,
                                   self.family,
                                   self.template,
                                   self.estimator,
                                   self.reverse,
                                   self.is_forward,
                                   self.electrons
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
        if self.specific_collider is not None: string += 'specific_collider={0!r}, '.format(self.specific_collider)
        if self.kinetics is not None: string += 'kinetics={0!r}, '.format(self.kinetics)
        if not self.reversible: string += 'reversible={0}, '.format(self.reversible)
        if self.transition_state is not None: string += 'transition_state={0!r}, '.format(self.transition_state)
        if self.duplicate: string += 'duplicate={0}, '.format(self.duplicate)
        if self.degeneracy != 1: string += 'degeneracy={0:.1f}, '.format(self.degeneracy)
        if self.pairs is not None: string += 'pairs={0}, '.format(self.pairs)
        if self.family: string += "family='{}', ".format(self.family)
        if self.template: string += "template={}, ".format(self.template)
        if self.electrons: string += "electrons={}, ".format(self.electrons)
        if self.comment != '': string += 'comment={0!r}, '.format(self.comment)
        string = string[:-2] + ')'
        return string

    def get_source(self):
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
        other.specific_collider = self.specific_collider
        other.degeneracy = self.degeneracy
        other.kinetics = deepcopy(self.kinetics)
        other.reversible = self.reversible
        other.transition_state = deepcopy(self.transition_state)
        other.duplicate = self.duplicate
        other.pairs = deepcopy(self.pairs)
        other.electrons = self.electrons

        # added for TemplateReaction information
        other.family = self.family
        other.template = self.template
        other.estimator = self.estimator
        other.reverse = self.reverse
        other.is_forward = self.is_forward

        return other

    def apply_solvent_correction(self, solvent):
        """
        apply kinetic solvent correction in this case the parameters are dGTSsite instead of GTS
        """
        from rmgpy.data.rmg import get_db
        solvation_database = get_db('solvation')
        solvent_data = solvation_database.get_solvent_data(solvent)
        
        
        if isinstance(self.kinetics, Marcus):
            solvent_struct = solvation_database.get_solvent_structure(solvent)[0]
            solv_solute_data = solvation_database.get_solute_data(solvent_struct.copy(deep=True))
            Rsolv = math.pow((75 * solv_solute_data.V / constants.pi / constants.Na),
                          (1.0 / 3.0)) / 100
            Rtot = 0.0
            Ner = 0
            Nep = 0
            for spc in self.reactants:
                spc_solute_data = solvation_database.get_solute_data(spc.copy(deep=True))
                spc_solute_data.set_mcgowan_volume(spc)
                R = math.pow((75 * spc_solute_data.V / constants.pi / constants.Na),
                            (1.0 / 3.0)) / 100
                Rtot += R 
                Ner += spc.get_net_charge()
            for spc in self.products:
                Nep += spc.get_net_charge()
            
            Rtot += Rsolv #radius of reactants plus first solvation shell
            self.lmbd_o = constants.Na*(constants.e*(Nep-Ner))**2/(8.0*constants.pi*constants.epsilon_0*Rtot)*(1.0/solvent_data.n**2 - 1.0/solvent_data.eps)
            return
        
        site_data = to_soluteTSdata(self.kinetics.solute)

        #compute x from gas phase
        GR = 0.0
        GP = 0.0
        for reactant in self.reactants:
            try:
                GR += reactant.get_free_energy(298.0)
            except Exception:
                logging.error("Problem with reactant {!r} in reaction {!s}".format(reactant, self))
                raise
        for product in self.products:
            try:
                GP += product.get_free_energy(298.0)
            except Exception:
                logging.error("Problem with product {!r} in reaction {!s}".format(reactant, self))
                raise
        
        GTS = self.kinetics.Ea.value_si + GR

        #x = abs(GTS - GR) / (abs(GP - GTS) + abs(GR - GTS))
        dGrxn = GP-GR
        if dGrxn > 0:
            x = 1.0
        else:
            x = 0.0

        dHR = 0.0
        dSR = 0.0
        for spc in self.reactants:
            spc_solute_data = solvation_database.get_solute_data(spc.copy(deep=True))
            spc_soluteTS_data = to_soluteTSdata(spc_solute_data)
            site_data += spc_soluteTS_data*(1.0-x)
            spc_correction = solvation_database.get_solvation_correction(spc_solute_data, solvent_data)
            dHR += spc_correction.enthalpy
            dSR += spc_correction.entropy

        for spc in self.products:
            spc_solute_data = to_soluteTSdata(solvation_database.get_solute_data(spc.copy(deep=True)))
            site_data += spc_solute_data*x

        dGTS,dHTS = site_data.calculate_corrections(solvent_data)
        dSTS = (dHTS - dGTS)/298.0

        dH = dHTS-dHR
        dA = np.exp((dSTS-dSR)/constants.R)
        self.kinetics.Ea.value_si += dH
        self.kinetics.A.value_si *= dA
        self.kinetics.comment += "\nsolvation correction raised barrier by {0} kcal/mol and prefactor by factor of {1}".format(dH/4184.0,dA)

################################################################################

class ReactionRecipe(object):
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

    def add_action(self, action):
        """
        Add an `action` to the reaction recipe, where `action` is a list
        containing the action name and the required parameters, as indicated in
        the table above.
        """
        self.actions.append(action)

    def get_reverse(self):
        """
        Generate a reaction recipe that, when applied, does the opposite of
        what the current recipe does, i.e., it is the recipe for the reverse
        of the reaction that this is the recipe for.
        """
        other = ReactionRecipe()
        for action in reversed(self.actions):  # Play the reverse recipe in the reverse order
            if action[0] == 'CHANGE_BOND':
                other.add_action(['CHANGE_BOND', action[1], str(-int(action[2])), action[3]])
            elif action[0] == 'FORM_BOND':
                other.add_action(['BREAK_BOND', action[1], action[2], action[3]])
            elif action[0] == 'BREAK_BOND':
                other.add_action(['FORM_BOND', action[1], action[2], action[3]])
            elif action[0] == 'LOSE_RADICAL':
                other.add_action(['GAIN_RADICAL', action[1], action[2]])
            elif action[0] == 'GAIN_RADICAL':
                other.add_action(['LOSE_RADICAL', action[1], action[2]])
            elif action[0] == 'GAIN_CHARGE':
                other.add_action(['LOSE_CHARGE', action[1], action[2]])
            elif action[0] == 'LOSE_CHARGE':
                other.add_action(['GAIN_CHARGE', action[1], action[2]])
            elif action[0] == 'LOSE_PAIR':
                other.add_action(['GAIN_PAIR', action[1], action[2]])
            elif action[0] == 'GAIN_PAIR':
                other.add_action(['LOSE_PAIR', action[1], action[2]])
        return other

    def _apply(self, struct, forward, unique):
        """
        Apply the reaction recipe to the set of molecules contained in
        `structure`, a single Structure object that contains one or more
        structures. The `forward` parameter is used to indicate
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
                struct.reset_connectivity_values()

                label1, info, label2 = action[1:]

                if label1 != label2:
                    # Find associated atoms
                    atom1 = struct.get_labeled_atoms(label1)[0]
                    atom2 = struct.get_labeled_atoms(label2)[0]
                else:
                    atoms = struct.get_labeled_atoms(label1)  # should never have more than two if this action is valid
                    if len(atoms) > 2:
                        raise InvalidActionError('Invalid atom labels encountered.')
                    atom1, atom2 = atoms

                if atom1 is None or atom2 is None or atom1 is atom2:
                    raise InvalidActionError('Invalid atom labels encountered.')

                # Apply the action
                if action[0] == 'CHANGE_BOND':
                    info = int(info)
                    # Check first to see if we have a bond
                    if not struct.has_bond(atom1, atom2):
                        if info < 1:
                            raise InvalidActionError('Attempted to change a nonexistent bond.')
                        # If we do not have a bond, it might be because we are trying to change a vdW bond
                        # Lets see if one of that atoms is a surface site,
                        # If we have a surface site, we will make a single bond, then change it by info - 1
                        is_vdW_bond = False
                        for atom in (atom1, atom2):
                            if atom.is_surface_site():
                                is_vdW_bond = True
                                break
                        if not is_vdW_bond: # no surface site, so no vdW bond
                            raise InvalidActionError('Attempted to change a nonexistent bond.')
                        else: # we found a surface site, so we will make a single bond
                            bond = GroupBond(atom1, atom2, order=[1]) if pattern else Bond(atom1, atom2, order=1)
                            struct.add_bond(bond)
                            atom1.apply_action(['FORM_BOND', label1, 1, label2])
                            atom2.apply_action(['FORM_BOND', label1, 1, label2])
                            # Now subtract 1 from info
                            info -= 1
                            # If info is 0, then we can continue since we don't have to change the bond
                            if info == 0:
                                continue
                    bond = struct.get_bond(atom1, atom2)
                    if bond.is_benzene():
                        struct.props['validAromatic'] = False
                    if forward:
                        atom1.apply_action(['CHANGE_BOND', label1, info, label2])
                        atom2.apply_action(['CHANGE_BOND', label1, info, label2])
                        bond.apply_action(['CHANGE_BOND', label1, info, label2])
                        if pattern:
                            if bond.is_van_der_waals():
                                if atom1.is_surface_site():
                                    atom1.atomtype = [ATOMTYPES['Xv']]
                                else:
                                    atom2.atomtype = [ATOMTYPES['Xv']]
                    else:
                        atom1.apply_action(['CHANGE_BOND', label1, -info, label2])
                        atom2.apply_action(['CHANGE_BOND', label1, -info, label2])
                        bond.apply_action(['CHANGE_BOND', label1, -info, label2])
                        if pattern:
                            if bond.is_van_der_waals():
                                if atom1.is_surface_site():
                                    atom1.atomtype = [ATOMTYPES['Xv']]
                                else:
                                    atom2.atomtype = [ATOMTYPES['Xv']]
                elif (action[0] == 'FORM_BOND' and forward) or (action[0] == 'BREAK_BOND' and not forward):
                    # Form bond between atom1 and atom2
                    if struct.has_bond(atom1, atom2):
                        raise InvalidActionError('Attempted to create an existing bond.')
                    if info not in (1, 0):  # Can only form single or vdW bonds
                        raise InvalidActionError('Attempted to create bond of type {:!r}'.format(info))
                    # we need to make sure we are not forming
                    # a surface bond to an atom that is already bonded to the surface
                    if atom1.is_surface_site() and atom2.is_bonded_to_surface():
                        raise InvalidActionError('Attempted to form a surface bond to an atom already bonded to surface.')
                    elif atom2.is_surface_site() and atom1.is_bonded_to_surface():
                        raise InvalidActionError('Attempted to form a surface bond to an atom already bonded to surface.')
                    bond = GroupBond(atom1, atom2, order=[info]) if pattern else Bond(atom1, atom2, order=info)
                    struct.add_bond(bond)
                    atom1.apply_action(['FORM_BOND', label1, info, label2])
                    atom2.apply_action(['FORM_BOND', label1, info, label2])
                elif (action[0] == 'BREAK_BOND' and forward) or (action[0] == 'FORM_BOND' and not forward):
                    # Break bond between atom1 and atom2
                    if not struct.has_bond(atom1, atom2):
                        if info == 0:
                            if atom1.is_surface_site() or atom2.is_surface_site():
                                # We are trying to break a vdW bond, but the atoms are not connected in
                                # the graph. The bond will break when we split the merged products
                                # in the `apply_recipe()` functions. Thus, there is nothing to do here,
                                # so we continue to the next action.
                                continue
                        raise InvalidActionError('Attempted to remove a nonexistent bond.')
                    bond = struct.get_bond(atom1, atom2)
                    struct.remove_bond(bond)
                    atom1.apply_action(['BREAK_BOND', label1, info, label2])
                    atom2.apply_action(['BREAK_BOND', label1, info, label2])

            elif action[0] in ['LOSE_RADICAL', 'GAIN_RADICAL', 'LOSE_CHARGE', 'GAIN_CHARGE']:

                label, change = action[1:]
                change = int(change)

                # Find associated atom
                atoms = struct.get_labeled_atoms(label)
                for atom in atoms:
                    if atom is None:
                        raise InvalidActionError('Unable to find atom with label "{0}" while applying '
                                                 'reaction recipe.'.format(label))

                    # Apply the action
                    for i in range(change):
                        if (action[0] == 'GAIN_RADICAL' and forward) or (action[0] == 'LOSE_RADICAL' and not forward):
                            atom.apply_action(['GAIN_RADICAL', label, 1])
                        elif (action[0] == 'LOSE_RADICAL' and forward) or (action[0] == 'GAIN_RADICAL' and not forward):
                            atom.apply_action(['LOSE_RADICAL', label, 1])
                        elif (action[0] == 'LOSE_CHARGE' and forward) or (action[0] == 'GAIN_CHARGE' and not forward):
                            atom.apply_action(['LOSE_CHARGE', label, 1])
                        elif (action[0] == 'GAIN_CHARGE' and forward) or (action[0] == 'LOSE_CHARGE' and not forward):
                            atom.apply_action(['GAIN_CHARGE', label, 1])

            elif action[0] in ['LOSE_PAIR', 'GAIN_PAIR']:

                label, change = action[1:]
                change = int(change)

                # Find associated atom
                atoms = struct.get_labeled_atoms(label)

                for atom in atoms:
                    if atom is None:
                        raise InvalidActionError('Unable to find atom with label "{0}" while applying '
                                                 'reaction recipe.'.format(label))

                    # Apply the action
                    for i in range(change):
                        if (action[0] == 'GAIN_PAIR' and forward) or (action[0] == 'LOSE_PAIR' and not forward):
                            atom.apply_action(['GAIN_PAIR', label, 1])
                        elif (action[0] == 'LOSE_PAIR' and forward) or (action[0] == 'GAIN_PAIR' and not forward):
                            atom.apply_action(['LOSE_PAIR', label, 1])

            else:
                raise InvalidActionError('Unknown action "' + action[0] + '" encountered.')

    def apply_forward(self, struct, unique=True):
        """
        Apply the forward reaction recipe to `molecule`, a single
        :class:`Molecule` object.
        """
        return self._apply(struct, True, unique)

    def apply_reverse(self, struct, unique=True):
        """
        Apply the reverse reaction recipe to `molecule`, a single
        :class:`Molecule` object.
        """
        return self._apply(struct, False, unique)


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
    `reversible`        ``Boolean``                     Is family reversible? (True by default)
    `forward_template`  :class:`Reaction`               The forward reaction template
    `forward_recipe`    :class:`ReactionRecipe`         The steps to take when applying the forward reaction to a set of reactants
    `reverse_template`  :class:`Reaction`               The reverse reaction template
    `reverse_recipe`    :class:`ReactionRecipe`         The steps to take when applying the reverse reaction to a set of reactants
    `forbidden`         :class:`ForbiddenStructures`    (Optional) Forbidden product structures in either direction
    `own_reverse`       ``Boolean``                     It's its own reverse?
    'boundary_atoms'    list                            Labels which define the boundaries of end groups in backbone/end families
    `tree_distances`    dict                            The default distance from parent along each tree, if not set default is 1 for every tree
    'save_order'        ``Boolean``                     Whether to preserve atom order when manipulating structures.
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
                 short_desc='',
                 long_desc='',
                 forward_template=None,
                 forward_recipe=None,
                 reverse_template=None,
                 reverse_recipe=None,
                 forbidden=None,
                 boundary_atoms=None,
                 tree_distances=None,
                 save_order=False,
                 ):
        Database.__init__(self, entries, top, label, name, short_desc, long_desc)
        self.reverse = reverse
        self.reversible = reversible
        self.forward_template = forward_template
        self.forward_recipe = forward_recipe
        self.reverse_template = reverse_template
        self.reverse_recipe = reverse_recipe
        self.forbidden = forbidden
        self.own_reverse = forward_template is not None and reverse_template is None
        self.boundary_atoms = boundary_atoms
        self.tree_distances = tree_distances
        self.save_order = save_order

        # Kinetics depositories of training and test data
        self.groups = None
        self.rules = None
        self.depositories = []

    def __repr__(self):
        return '<ReactionFamily "{0}">'.format(self.label)

    def distribute_tree_distances(self):
        """
        Fills in nodal_distance (the distance between an entry and its parent)
        if not already entered with the value from tree_distances associated
        with the tree the entry comes from
        """
        tree_distances = self.tree_distances
        top_labels = [i.label for i in self.groups.top]

        if len(top_labels) != len(tree_distances):
            raise ValueError('tree_distances does not have the same number of '
                             'entries as there are top nodes in the family')

        for entry_name, entry in self.groups.entries.items():
            top_entry = entry
            while not (top_entry.parent is None):  # get the top for the tree entry is in
                top_entry = top_entry.parent
            if top_entry.label in top_labels:  # filtering out product nodes
                if entry.nodal_distance is None:
                    entry.nodal_distance = tree_distances[top_entry.label]

    def load(self, path, local_context=None, global_context=None, depository_labels=None):
        """
        Load a kinetics database from a file located at `path` on disk.

        If `depository_labels` is a list, eg. ['training','PrIMe'], then only those
        depositories are loaded, and they are searched in that order when
        generating kinetics.

        If depository_labels is None then load 'training' first then everything else.
        If depository_labels is not None then load in the order specified in depository_labels.
        """
        local_context['recipe'] = self.load_recipe
        local_context['template'] = self.load_template
        local_context['forbidden'] = self.load_forbidden
        local_context['True'] = True
        local_context['False'] = False
        local_context['reverse'] = None
        local_context['reversible'] = None
        local_context['boundaryAtoms'] = None
        local_context['treeDistances'] = None
        local_context['reverseMap'] = None
        local_context['reactantNum'] = None
        local_context['productNum'] = None
        local_context['autoGenerated'] = False
        local_context['allowChargedSpecies'] = False
        local_context['electrons'] = 0
        self.groups = KineticsGroups(label='{0}/groups'.format(self.label))
        logging.debug("Loading kinetics family groups from {0}".format(os.path.join(path, 'groups.py')))
        Database.load(self.groups, os.path.join(path, 'groups.py'), local_context, global_context)
        self.name = self.label
        self.boundary_atoms = local_context.get('boundaryAtoms', None)
        self.tree_distances = local_context.get('treeDistances', None)
        self.reverse_map = local_context.get('reverseMap', None)

        self.reactant_num = local_context.get('reactantNum', None)
        self.product_num = local_context.get('productNum', None)

        self.auto_generated = local_context.get('autoGenerated', False)
        self.allow_charged_species = local_context.get('allowChargedSpecies', False)
        self.electrons = local_context.get('electrons', 0)

        if self.reactant_num:
            self.groups.reactant_num = self.reactant_num
        else:
            self.groups.reactant_num = len(self.forward_template.reactants)

        # Generate the reverse template if necessary
        self.forward_template.reactants = [self.groups.entries[label] for label in self.forward_template.reactants]
        if self.own_reverse:
            self.forward_template.products = self.forward_template.reactants[:]
            self.reverse_template = None
            self.reverse_recipe = self.forward_recipe.get_reverse()
        else:
            self.reverse = local_context.get('reverse', None)
            self.reversible = True if local_context.get('reversible', None) is None else local_context.get('reversible', None)
            self.forward_template.products = self.generate_product_template(self.forward_template.reactants)
            for entry in self.forward_template.products:
                if isinstance(entry.item,Group):
                    entry.item.update()
            if self.reversible:
                self.reverse_template = Reaction(reactants=self.forward_template.products,
                                                 products=self.forward_template.reactants)
                self.reverse_recipe = self.forward_recipe.get_reverse()
                if self.reverse is None:
                    self.reverse = '{0}_reverse'.format(self.label)

        self.rules = KineticsRules(label='{0}/rules'.format(self.label),auto_generated=self.auto_generated)
        logging.debug("Loading kinetics family rules from {0}".format(os.path.join(path, 'rules.py')))
        self.rules.load(os.path.join(path, 'rules.py'), local_context, global_context)

        # load the groups indicated in the entry label
        for label, entries in self.rules.entries.items():
            nodes = label.split(';')
            reactants = [self.groups.entries[node] for node in nodes]
            reaction = Reaction(reactants=reactants, products=[])
            for entry in entries:
                entry.item = reaction
        self.depositories = []

        top_labels = [i.label for i in self.groups.top]
        if self.tree_distances is None:
            self.tree_distances = {top_entry: 1 for top_entry in top_labels}

        self.distribute_tree_distances()

        if depository_labels == 'all':
            # Load everything. This option is generally used for working with the database
            # load all the remaining depositories, in order returned by os.walk
            for root, dirs, files in os.walk(path):
                for name in dirs:
                    # if not f.endswith('.py'): continue
                    # name = f.split('.py')[0]
                    # if name not in ['groups', 'rules']:
                    f_path = os.path.join(path, name, 'reactions.py')
                    label = '{0}/{1}'.format(self.label, name)
                    depository = KineticsDepository(label=label)
                    logging.debug("Loading kinetics family depository from {0}".format(f_path))
                    depository.load(f_path, local_context, global_context)
                    self.depositories.append(depository)

            return

        if not depository_labels:
            # If depository labels is None or there are no depositories listed, then use the training
            # depository and add them to the RMG rate rules by default:
            depository_labels = ['training']
        if depository_labels:
            # If there are depository labels, load them in the order specified, but
            # append the training reactions unless the user specifically declares it not
            # to be included with a '!training' flag
            if '!training' not in depository_labels:
                if 'training' not in depository_labels:
                    depository_labels.append('training')

        for name in depository_labels:
            if name == '!training':
                continue
            label = '{0}/{1}'.format(self.label, name)
            # f = name+'.py'
            f_path = os.path.join(path, name, 'reactions.py')
            if not os.path.exists(f_path):
                logging.warning("Requested depository {0} does not exist".format(f_path))
                continue
            depository = KineticsDepository(label=label)
            logging.debug("Loading kinetics family depository from {0}".format(f_path))
            depository.load(f_path, local_context, global_context)
            self.depositories.append(depository)

    def load_template(self, reactants, products, ownReverse=False):
        """
        Load information about the reaction template.
        Note that argument names are retained for backward compatibility
        with loading database files.
        """
        self.forward_template = Reaction(reactants=reactants, products=products)
        self.own_reverse = ownReverse

    def load_recipe(self, actions):
        """
        Load information about the reaction recipe.
        """
        # Remaining lines are reaction recipe for forward reaction
        self.forward_recipe = ReactionRecipe()
        for action in actions:
            action[0] = action[0].upper()
            valid_actions = [
                'CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND', 'GAIN_RADICAL', 'LOSE_RADICAL',
                'GAIN_CHARGE', 'LOSE_CHARGE', 'GAIN_PAIR', 'LOSE_PAIR'
            ]
            if action[0] not in valid_actions:
                raise InvalidActionError('Action {0} is not a recognized action. '
                                         'Should be one of {1}'.format(actions[0], valid_actions))
            self.forward_recipe.add_action(action)

    def load_forbidden(self, label, group, shortDesc='', longDesc=''):
        """
        Load information about a forbidden structure.
        Note that argument names are retained for backward compatibility
        with loading database files.
        """
        if not self.forbidden:
            self.forbidden = ForbiddenStructures()
        self.forbidden.load_entry(label=label, group=group, shortDesc=shortDesc, longDesc=longDesc)

    def save_entry(self, f, entry):
        """
        Write the given `entry` in the thermo database to the file object `f`.
        """
        return save_entry(f, entry)

    def save_training_reactions(self, reactions, reference=None, reference_type='', short_desc='', long_desc='',
                                rank=3):
        """
        This function takes a list of reactions appends it to the training reactions file.  It ignores the existence of
        duplicate reactions.

        The rank for each new reaction's kinetics is set to a default value of 3 unless the user specifies differently
        for those reactions.

        For each entry, the long description is imported from the kinetics comment.
        """

        if not isinstance(reference, list):
            reference = [reference] * len(reactions)
        if not isinstance(reference_type, list):
            reference_type = [reference_type] * len(reactions)
        if not isinstance(short_desc, list):
            short_desc = [short_desc] * len(reactions)
        if not isinstance(long_desc, list):
            long_desc = [long_desc] * len(reactions)
        if not isinstance(rank, list):
            rank = [rank] * len(reactions)

        training_path = os.path.join(settings['database.directory'], 'kinetics', 'families',
                                     self.label, 'training')

        dictionary_path = os.path.join(training_path, 'dictionary.txt')

        # Load the old set of the species of the training reactions
        species_dict = Database().get_species(dictionary_path)

        # Add new unique species with labeled atoms into species_dict
        for rxn in reactions:
            for spec in (rxn.reactants + rxn.products):
                for ex_spec in species_dict.values():
                    if ex_spec.molecule[0].get_formula() != spec.molecule[0].get_formula():
                        continue
                    else:
                        spec_labeled_atoms = spec.molecule[0].get_all_labeled_atoms()
                        spcs_labels = sorted(spec_labeled_atoms.keys())
                        ex_spec_labeled_atoms = ex_spec.molecule[0].get_all_labeled_atoms()
                        ex_spcs_labels = sorted(ex_spec_labeled_atoms.keys())
                        if spcs_labels != ex_spcs_labels:
                            # the species have different labels, therefore not a match
                            continue
                        initial_map = {}
                        try:
                            for atomLabel in spec_labeled_atoms:
                                initial_map[spec_labeled_atoms[atomLabel]] = ex_spec_labeled_atoms[atomLabel]
                        except KeyError:
                            # Atom labels did not match, therefore not a match
                            continue
                        if spec.molecule[0].is_isomorphic(ex_spec.molecule[0], initial_map):
                            spec.label = ex_spec.label
                            break
                else:  # No isomorphic existing species found
                    spec_formula = spec.molecule[0].get_formula()
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
            depository = self.get_training_depository()
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
            index = max_index + i + 1
            entry = Entry(
                index=index,
                label=str(reaction),
                item=reaction,
                data=reaction.kinetics,
                reference=reference[i],
                reference_type=reference_type[i],
                short_desc=str(short_desc[i]),
                long_desc=str(long_desc[i]),
                rank=rank[i],
            )

            # Add this entry to the loaded depository so it is immediately usable
            depository.entries[index] = entry
            # Write the entry to the reactions.py file
            self.save_entry(training_file, entry)

        training_file.close()

        # save species to dictionary
        with open(dictionary_path, 'w') as f:
            for label in species_dict.keys():
                f.write(species_dict[label].molecule[0].to_adjacency_list(label=label, remove_h=False))
                f.write('\n')

    def save(self, path):
        """
        Save the current database to the file at location `path` on disk.
        """
        self.save_groups(os.path.join(path, 'groups.py'))
        self.rules.save(os.path.join(path, 'rules.py'))
        for depository in self.depositories:
            self.save_depository(depository, os.path.join(path, '{0}'.format(depository.label[len(self.label) + 1:])))

    def save_depository(self, depository, path):
        """
        Save the given kinetics family `depository` to the location `path` on
        disk.
        """
        depository.save_dictionary(os.path.join(path, 'dictionary.txt'))
        depository.save(os.path.join(path, 'reactions.py'))

    def save_groups(self, path):
        """
        Save the current database to the file at location `path` on disk.
        """
        entries = self.groups.get_entries_to_save()

        # Write the header
        f = codecs.open(path, 'w', 'utf-8')
        f.write('#!/usr/bin/env python\n')
        f.write('# encoding: utf-8\n\n')
        f.write('name = "{0}/groups"\n'.format(self.name))
        f.write('shortDesc = "{0}"\n'.format(self.groups.short_desc))
        f.write('longDesc = """\n')
        f.write(self.groups.long_desc)
        f.write('\n"""\n\n')

        # Write the template
        f.write('template(reactants=[{0}], products=[{1}], ownReverse={2})\n\n'.format(
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forward_template.reactants]),
            ', '.join(['"{0}"'.format(entry.label) for entry in self.forward_template.products]),
            self.own_reverse))

        # Write reverse name
        if not self.own_reverse:
            if self.reverse is not None:
                f.write('reverse = "{0}"\n'.format(self.reverse))
            else:
                f.write('reverse = None\n')

        f.write('reversible = {0}\n\n'.format(self.reversible))

        if self.reverse_map is not None:
            f.write('reverseMap = {0}\n\n'.format(self.reverse_map))

        if self.reactant_num is not None:
            f.write('reactantNum = {0}\n\n'.format(self.reactant_num))
        if self.product_num is not None:
            f.write('productNum = {0}\n\n'.format(self.product_num))

        if self.auto_generated is not None:
            f.write('autoGenerated = {0}\n\n'.format(self.auto_generated))

        if self.allow_charged_species:
            f.write('allowChargedSpecies = {0}\n\n'.format(self.allow_charged_species))

        if self.electrons != 0:
            f.write('electrons = {0}\n\n'.format(self.electrons))

        # Write the recipe
        f.write('recipe(actions=[\n')
        for action in self.forward_recipe.actions:
            f.write('    {0!r},\n'.format(action))
        f.write('])\n\n')

        if self.boundary_atoms:
            f.write('boundaryAtoms = ["{0}", "{1}"]'.format(self.boundary_atoms[0], self.boundary_atoms[1]))
            f.write('\n\n')

        # Save the entries
        for entry in entries:
            self.save_entry(f, entry)

        # Write the tree
        if len(self.groups.top) > 0:
            f.write('tree(\n')
            f.write('"""\n')
            f.write(self.generate_old_tree(self.groups.top, 1))
            f.write('"""\n')
            f.write(')\n\n')

        # Save forbidden structures, if present
        if self.forbidden is not None:
            entries = list(self.forbidden.entries.values())
            entries.sort(key=lambda x: x.label)
            for entry in entries:
                self.forbidden.save_entry(f, entry, name='forbidden')

        f.close()

    def generate_product_template(self, reactants0):
        """
        Generate the product structures by applying the reaction template to
        the top-level nodes. For reactants defined by multiple structures, only
        the first is used here; it is assumed to be the most generic.
        """

        # First, generate a list of reactant structures that are actual
        # structures, rather than unions
        reactant_structures = []
        for reactant in reactants0:
            if isinstance(reactant, list):
                reactants = [reactant[0]]
            else:
                reactants = [reactant]

            for s in reactants:
                struct = s.item
                if isinstance(struct, LogicNode):
                    all_structures = struct.get_possible_structures(self.groups.entries)
                    reactant_structures.append(all_structures)
                else:
                    reactant_structures.append([struct])

        # Second, get all possible combinations of reactant structures
        reactant_structures = get_all_combinations(reactant_structures)

        # Third, generate all possible product structures by applying the
        # recipe to each combination of reactant structures
        # Note that bimolecular products are split by labeled atoms
        product_structures = []
        for reactant_structure in reactant_structures:
            product_structure = self.apply_recipe(reactant_structure, forward=True, unique=False)
            if product_structure:
                product_structures.append(product_structure)

        # Fourth, remove duplicates from the lists
        product_structure_list = [[] for i in range(len(product_structures[0]))]
        for product_structure in product_structures:
            for i, struct in enumerate(product_structure):
                for s in product_structure_list[i]:
                    try:
                        if s.is_identical(struct): break
                    except KeyError:
                        logging.error(struct.to_adjacency_list())
                        logging.error(s.to_adjacency_list())
                        raise
                else:
                    product_structure_list[i].append(struct)
        # Fifth, associate structures with product template
        product_set = []
        for index, products in enumerate(product_structure_list):
            label = self.forward_template.products[index]
            if len(products) == 1:
                entry = Entry(
                    label=label,
                    item=products[0],
                )
                self.groups.entries[entry.label] = entry
                product_set.append(entry)
            else:
                children = []
                counter = 0
                for product in products:
                    entry = Entry(
                        label='{0}{1:d}'.format(label, counter + 1),
                        item=product,
                    )
                    children.append(entry)
                    self.groups.entries[entry.label] = entry
                    counter += 1

                # Enter the parent of the groups as a logicOr of all the products
                entry = Entry(
                    label=label,
                    item=LogicOr([child.label for child in children], invert=False),
                    children=children,
                )
                self.groups.entries[entry.label] = entry
                # Make this entry the parent of all its children
                for child in children:
                    child.parent = entry
                counter += 1
                product_set.append(entry)

        return product_set

    def has_rate_rule(self, template):
        """
        Return ``True`` if a rate rule with the given `template` currently
        exists, or ``False`` otherwise.
        """
        return self.rules.has_rule(template)

    def get_rate_rule(self, template):
        """
        Return the rate rule with the given `template`. Raises a
        :class:`ValueError` if no corresponding entry exists.
        """
        entry = self.rules.get_rule(template)
        if entry is None:
            raise ValueError('No entry for template {0}.'.format(template))
        return entry

    def add_rules_from_training(self, thermo_database=None, train_indices=None):
        """
        For each reaction involving real reactants and products in the training
        set, add a rate rule for that reaction.
        """
        if self.auto_generated:
            warnings.warn(f'add_rules_from_training should be only called for non-ATG families, '
                          f'but {self.label} is an ATG family. Skip this function call. '
                          f'Calling add_rules_from_training on ATG families may be deprecated in RMG-Py 3.2',
                          DeprecationWarning)
            return

        try:
            depository = self.get_training_depository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Must be because you turned off the training depository.')
            return

        # Determine number of parallel processes.
        from rmgpy.rmg.main import determine_procnum_from_ram
        procnum = determine_procnum_from_ram()

        tentries = depository.entries

        index = max([e.index for e in self.rules.get_entries()] or [0]) + 1

        entries = list(depository.entries.values())
        entries.sort(key=lambda x: x.index)

        if train_indices is not None:
            entries = np.array(entries)
            entries = entries[train_indices]

        reverse_entries = []
        for entry in entries:
            try:
                template = self.get_reaction_template(entry.item)
            except UndeterminableKineticsError:
                # Some entries might be stored in the reverse direction for
                # this family; save them so we can try this
                reverse_entries.append(entry)
                continue

            tentries[entry.index].item.is_forward = True

            data = deepcopy(entry.data)
            data.change_t0(1)

            if type(data) is Arrhenius:
                # more specific than isinstance(data,Arrhenius) because we want to exclude inherited subclasses!
                data = data.to_arrhenius_ep()
            elif isinstance(data, StickingCoefficient):
                data = StickingCoefficientBEP(
                    # todo: perhaps make a method StickingCoefficient.StickingCoefficientBEP
                    #  analogous to Arrhenius.to_arrhenius_ep
                    A=deepcopy(data.A),
                    n=deepcopy(data.n),
                    alpha=0,
                    E0=deepcopy(data.Ea),
                    Tmin=deepcopy(data.Tmin),
                    Tmax=deepcopy(data.Tmax),
                    coverage_dependence=deepcopy(data.coverage_dependence),
                )
            elif isinstance(data, SurfaceArrhenius):
                data = SurfaceArrheniusBEP(
                    # todo: perhaps make a method SurfaceArrhenius.toSurfaceArrheniusBEP
                    #  analogous to Arrhenius.to_arrhenius_ep
                    A=deepcopy(data.A),
                    n=deepcopy(data.n),
                    alpha=0,
                    E0=deepcopy(data.Ea),
                    Tmin=deepcopy(data.Tmin),
                    Tmax=deepcopy(data.Tmax),
                    coverage_dependence=deepcopy(data.coverage_dependence),
                )
            elif isinstance(data, SurfaceChargeTransfer):
                for reactant in entry.item.reactants:
                    # Clear atom labels to avoid effects on thermo generation, ok because this is a deepcopy
                    reactant_copy = reactant.copy(deep=True)
                    reactant_copy.molecule[0].clear_labeled_atoms()
                    reactant_copy.generate_resonance_structures()
                    reactant.thermo = thermo_database.get_thermo_data(reactant_copy, training_set=True)
                for product in entry.item.products:
                    product_copy = product.copy(deep=True)
                    product_copy.molecule[0].clear_labeled_atoms()
                    product_copy.generate_resonance_structures()
                    product.thermo = thermo_database.get_thermo_data(product_copy, training_set=True)
                V = data.V0.value_si
                dGrxn = entry.item._get_free_energy_of_charge_transfer_reaction(298,V)
                data = data.to_surface_charge_transfer_bep(dGrxn,0.0)
            else:
                raise NotImplementedError("Unexpected training kinetics type {} for {}".format(type(data), entry))

            new_entry = Entry(
                index=index,
                label=';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template], products=[]),
                data=data,
                rank=entry.rank,
                reference=entry.reference,
                short_desc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.short_desc,
                long_desc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.long_desc,
            )
            new_entry.data.comment = "From training reaction {1} used for {0}".format(
                ';'.join([g.label for g in template]), entry.index)

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

            if not isinstance(entry.data, Arrhenius):
                print(self.label)
                assert False
            data = deepcopy(entry.data)
            data.change_t0(1)
            # Estimate the thermo for the reactants and products
            # training_set=True used later to does not allow species to match a liquid phase library
            # and get corrected thermo which will affect reverse rate calculation
            item = Reaction(reactants=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label)
                                       for m in entry.item.reactants],
                            products=[Species(molecule=[m.molecule[0].copy(deep=True)], label=m.label)
                                      for m in entry.item.products])

            if procnum > 1:
                # If QMTP and multiprocessing write QMTP files here in parallel.
                from rmgpy.rmg.input import get_input
                quantum_mechanics = get_input('quantum_mechanics')
                if quantum_mechanics:
                    quantum_mechanics.run_jobs(item.reactants + item.products, procnum=procnum)

            if entry.facet is None:
                metal = entry.metal # could be None
            else:
                metal = entry.metal + entry.facet

            for reactant in item.reactants:
                # Clear atom labels to avoid effects on thermo generation, ok because this is a deepcopy
                reactant.molecule[0].clear_labeled_atoms()
                reactant.generate_resonance_structures()
                reactant.thermo = thermo_database.get_thermo_data(reactant, training_set=True, metal_to_scale_to=metal)
            for product in item.products:
                product.molecule[0].clear_labeled_atoms()
                product.generate_resonance_structures()
                product.thermo = thermo_database.get_thermo_data(product, training_set=True, metal_to_scale_to=metal)
            # Now that we have the thermo, we can get the reverse k(T)
            item.kinetics = data
            data = item.generate_reverse_rate_coefficient()

            item = TemplateReaction(reactants=[m.molecule[0].copy(deep=True) for m in entry.item.products],
                                    products=[m.molecule[0].copy(deep=True) for m in entry.item.reactants])
            template = self.get_reaction_template(item)

            item.template = self.get_reaction_template_labels(item)
            new_degeneracy = self.calculate_degeneracy(item)

            if isinstance(entry.data, SurfaceArrhenius):
                data = SurfaceArrheniusBEP(
                    #  analogous to Arrhenius.to_arrhenius_ep
                    A=deepcopy(data.A),
                    n=deepcopy(data.n),
                    alpha=0,
                    E0=deepcopy(data.Ea),
                    Tmin=deepcopy(data.Tmin),
                    Tmax=deepcopy(data.Tmax),
                    coverage_dependence=deepcopy(data.coverage_dependence),
                )
            else:
                data = data.to_arrhenius_ep()

            new_entry = Entry(
                index=index,
                label=';'.join([g.label for g in template]),
                item=Reaction(reactants=[g.item for g in template],
                              products=[]),
                data=data,
                rank=entry.rank,
                reference=entry.reference,
                short_desc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.short_desc,
                long_desc="Rate rule generated from training reaction {0}. ".format(entry.index) + entry.long_desc,
                metal=entry.metal
            )
            new_entry.data.comment = "From training reaction {1} used for {0}".format(';'.join([g.label for g in template]), entry.index)

            new_entry.data.A.value_si /= new_degeneracy
            try:
                self.rules.entries[new_entry.label].append(new_entry)
            except KeyError:
                self.rules.entries[new_entry.label] = [new_entry]
            index += 1

    def get_root_template(self):
        """
        Return the root template for the reaction family. Most of the time this
        is the top-level nodes of the tree (as stored in the
        :class:`KineticsGroups` object), but there are a few exceptions (e.g.
        R_Recombination).
        """
        if len(self.forward_template.reactants) > len(self.groups.top):
            return self.forward_template.reactants
        else:
            return self.groups.top

    def fill_rules_by_averaging_up(self, verbose=False):
        """
        Fill in gaps in the kinetics rate rules by averaging child nodes
        recursively starting from the top level root template.
        """
        if self.auto_generated:
            warnings.warn(f'fill_rules_by_averaging_up should be only called for non-ATG families, '
                          f'but {self.label} is an ATG family. Skip this function call. '
                          f'Calling fill_rules_by_averaging_up on ATG families may be deprecated in RMG-Py 3.2',
                          DeprecationWarning)
            return
        self.rules.fill_rules_by_averaging_up(self.get_root_template(), {}, verbose)

    def apply_recipe(self, reactant_structures, forward=True, unique=True, relabel_atoms=True):
        """
        Apply the recipe for this reaction family to the list of
        :class:`Molecule` or :class:`Group` objects `reactant_structures`. The atoms
        of the reactant structures must already be tagged with the appropriate
        labels. Returns a list of structures corresponding to the products
        after checking that the correct number of products was produced.
        If ``relabel_atoms`` is ``True``, product atom labels of reversible families
        will be reversed to assist in identifying forbidden structures.
        """

        # There is some hardcoding of reaction families in this function, so
        # we need the label of the reaction family for this
        label = self.label.lower()

        # Merge reactant structures into single structure
        # Also copy structures so we don't modify the originals
        # Since the tagging has already occurred, both the reactants and the
        # products will have tags
        if any(isinstance(reactant, Fragment) for reactant in reactant_structures):
            reactant_structure = Fragment()
        elif isinstance(reactant_structures[0], Group):
            reactant_structure = Group()
        elif isinstance(reactant_structures[0], Molecule):
            reactant_structure = Molecule()
        for s in reactant_structures:
            reactant_structure = reactant_structure.merge(s.copy(deep=True))

        if forward:
            # Generate the product structure by applying the recipe
            self.forward_recipe.apply_forward(reactant_structure, unique)
        else:
            self.reverse_recipe.apply_forward(reactant_structure, unique)

        # Now that we have applied the recipe, let's start calling
        # this thing the product_structure (although it's the same object in memory)
        product_structure = reactant_structure

        if not product_structure.props['validAromatic']:
            if isinstance(product_structure, Molecule) or isinstance(product_structure, Fragment):
                # For molecules, kekulize the product to redistribute bonds appropriately
                product_structure.kekulize()
            else:
                # For groups, we ignore the product template for a purely aromatic group
                # If there is an analagous aliphatic group in the family, then the product template will be identical
                # There should NOT be any families that consist solely of aromatic reactant templates
                return []

        # If reaction family is its own reverse, relabel atoms
        # This allows comparison of the product species to forbidden
        #  structures which are labeled as reactants.
        # Unfortunately, this means that reaction family info is
        #  hardcoded, so this must be updated if the database changes.
        if not self.reverse_template and relabel_atoms:
            # Get atom labels for products
            atom_labels = {}
            for atom in product_structure.atoms:
                if atom.label != '':
                    atom_labels[atom.label] = atom

            if label in ('1,2_xy_interchange'):
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*4'
                atom_labels['*4'].label = '*1'

            if label in ('h_abstraction','f_abstraction','cl_abstraction','br_abstraction'):
                # '*2' is the H that migrates
                # it moves from '*1' to '*3'
                atom_labels['*1'].label = '*3'
                atom_labels['*3'].label = '*1'

            elif label == 'intra_h_migration':
                # '*3' is the H that migrates
                # swap the two ends between which the H moves
                atom_labels['*1'].label = '*2'
                atom_labels['*2'].label = '*1'
                # reverse all the atoms in the chain between *1 and *2
                highest = len(atom_labels)
                if highest > 4:
                    # swap *4 with *5
                    atom_labels['*4'].label = '*5'
                    atom_labels['*5'].label = '*4'
                if highest > 6:
                    # swap *6 with the highest, etc.
                    for i in range(6, highest + 1):
                        atom_labels['*{0:d}'.format(i)].label = '*{0:d}'.format(6 + highest - i)

            elif label == 'intra_ene_reaction':
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*2'
                atom_labels['*2'].label = '*1'
                atom_labels['*3'].label = '*5'
                atom_labels['*5'].label = '*3'

            elif label == '6_membered_central_c-c_shift':
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*3'
                atom_labels['*3'].label = '*1'
                atom_labels['*4'].label = '*6'
                atom_labels['*6'].label = '*4'

            elif label == '1,2_shiftc':
                # Labels for nodes are swapped
                atom_labels['*2'].label = '*3'
                atom_labels['*3'].label = '*2'

            elif label == 'intra_r_add_exo_scission':
                # Labels for nodes are swapped
                atom_labels['*1'].label = '*3'
                atom_labels['*3'].label = '*1'

            elif label == 'intra_substitutions_isomerization':
                # Swap *2 and *3
                atom_labels['*2'].label = '*3'
                atom_labels['*3'].label = '*2'

            elif label == 'surface_abstraction':
                atom_labels['*1'].label = '*3'
                atom_labels['*2'].label = '*5'
                atom_labels['*3'].label = '*1'
                atom_labels['*5'].label = '*2'

            elif label == 'surface_abstraction_single_vdw':
                # *3 migrates from *2-*1 to *4-*5
                # so swap *1 with *5, swap *2 with *4
                atom_labels['*1'].label = '*5'
                atom_labels['*5'].label = '*1'
                atom_labels['*2'].label = '*4'
                atom_labels['*4'].label = '*2'

            # Use the family's reverse map if it has one
            elif self.reverse_map:
                for label0,label1 in self.reverse_map.items():
                    atom_labels[label0] = label1

        if not forward:
            template = self.reverse_template
            product_num = self.reactant_num or len(template.products)
        else:
            template = self.forward_template
            product_num = self.product_num or len(template.products)

        # Split product structure into multiple species if necessary
        if self.auto_generated and isinstance(reactant_structures[0],Group) and self.product_num == 1:
            product_structures = [product_structure]
        else:
            product_structures = product_structure.split()

        # Make sure we've made the expected number of products
        if product_num != len(product_structures):
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

        # Remove vdW bonds
        for struct in product_structures:
            if isinstance(struct, Fragment):
                continue
            else:
                struct.remove_van_der_waals_bonds()

        # Make sure we don't create a different net charge between reactants and products
        reactant_net_charge = product_net_charge = 0
        for struc in reactant_structures:
            if isinstance(struc, Molecule):
                struc.update(sort_atoms=not self.save_order)
            else:
                struc.update()
            reactant_net_charge += struc.get_net_charge()


        is_molecule = True
        for struct in product_structures:
            # If product structures are Molecule objects, update their atom types
            # If product structures are Group objects and the reaction is in certain families
            # (families with charged substances), the charge of structures will be updated
            if isinstance(struct, Molecule):
                struct.update_charge()
                if isinstance(struct, Fragment):
                    struct.update()
                else:
                    struct.update(sort_atoms=not self.save_order)
            elif isinstance(struct, Group):
                is_molecule = False
                struct.reset_ring_membership()
                if label in ['1,2_insertion_co', 'r_addition_com', 'co_disproportionation',
                             'intra_no2_ono_conversion', 'lone_electron_pair_bond',
                             '1,2_nh3_elimination', '1,3_nh3_elimination']:
                    struct.update_charge()
            else:
                raise TypeError('Expecting Molecule or Group object, not {0}'.format(struct.__class__.__name__))
            product_net_charge += struct.get_net_charge()


        if self.electrons < 0:
            if forward:
                reactant_net_charge += self.electrons
            else:
                product_net_charge += self.electrons
        elif self.electrons > 0:
            if forward:
                product_net_charge -= self.electrons
            else:
                reactant_net_charge -= self.electrons

        if reactant_net_charge != product_net_charge and is_molecule:
            logging.debug(
                'The net charge of the reactants {0} differs from the net charge of the products {1} in reaction '
                'family {2}. Not generating this reaction.'.format(reactant_net_charge, product_net_charge, self.label))
            return None

        # If there are two product structures, place the one containing '*1' first
        if len(product_structures) == 2:
            if not product_structures[0].contains_labeled_atom('*1') and \
                    product_structures[1].contains_labeled_atom('*1'):
                product_structures.reverse()
        # If there are three product structures, sort them based on the lowest number label in each structure
        elif len(product_structures) == 3:
            lowest_labels = []
            for struct in product_structures:
                # Extract digits from labels and convert others (e.g., "*") to empty strings
                labels = [''.join(c for c in label if c.isdigit()) for label in struct.get_all_labeled_atoms().keys()]
                # Convert digits to integers and remove empty strings
                labels = [int(label) for label in labels if label]
                lowest_labels.append(min(labels))
            product_structures = [s for _, s in sorted(zip(lowest_labels, product_structures))]

        # If the template restricts multiplicity, we need to make sure that
        # a product matches the multiplicity-constrained template
        # because the template does not ensure that the
        # multiplicity restriction is obeyed
        if isinstance(reactant_structure, Molecule):
            if forward:
                template_groups = self.forward_template.products
            else:
                template_groups = self.forward_template.reactants
            for template in template_groups:
                # iterate through the template reactants and check to see if they have a multiplicity constraint
                if isinstance(template.item, Group):
                    if template.item.multiplicity != []:
                        # this template restricts multiplicity and needs to be checked
                        for struct in product_structures:
                            # iterate through the product structures to make sure that
                            # it matches a mulitplicity-constrained template reactant
                            match = self._match_reactant_to_template(struct, template)
                            if match:
                                # A product structure matches the template!
                                break
                        if not match:
                            # No product matched the template reactant
                            # Therefore, this reaction is invalid
                            if logging.getLogger().isEnabledFor(logging.DEBUG):
                                logging.debug(
                                    'No product structures matched %s which has a multiplicity of %r\n'
                                    'Template:\n%s\n'
                                    'Product structures:\n',
                                    template, template.item.multiplicity, template.item.to_adjacency_list() )
                                for struct in product_structures:
                                    logging.debug(f'{struct}\n{struct.to_adjacency_list()}\n')
                            return None

        # Return the product structures
        return product_structures

    def _generate_product_structures(self, reactant_structures, maps, forward, relabel_atoms=True):
        """
        For a given set of `reactant_structures` and a given set of `maps`,
        generate and return the corresponding product structures. The
        `reactant_structures` parameter should be given in the order the
        reactants are stored in the reaction family template. The `maps`
        parameter is a list of mappings of the top-level tree node of each
        *template* reactant to the corresponding *structure*. This function
        returns a list of the product structures.
        If ``relabel_atoms`` is ``True``, product atom labels of reversible families
        will be reversed to assist in identifying forbidden structures.
        """

        # Clear any previous atom labeling from all reactant structures
        for struct in reactant_structures:
            struct.clear_labeled_atoms()

        # Tag atoms with labels
        for m in maps:
            for reactant_atom, template_atom in m.items():
                reactant_atom.label = template_atom.label

        # Check that reactant structures are allowed in this family
        # If not, then stop
        for struct in reactant_structures:
            if self.is_molecule_forbidden(struct):
                raise ForbiddenStructureException()

        # Generate the product structures by applying the forward reaction recipe
        try:
            product_structures = self.apply_recipe(reactant_structures, forward=forward, relabel_atoms=relabel_atoms)
            if not product_structures:
                return None
        except (InvalidActionError, KekulizationError, AtomTypeError):
            # If unable to apply the reaction recipe, then return no product structures
            return None
        except ActionError:
            logging.error('Could not generate product structures for reaction family {0} in {1} '
                          'direction'.format(self.label, 'forward' if forward else 'reverse'))
            logging.info('Reactant structures:')
            for struct in reactant_structures:
                logging.info('{0}\n{1}\n'.format(struct, struct.to_adjacency_list()))
            raise

        # Apply the generated species constraints (if given)
        for struct in product_structures:
            if self.is_molecule_forbidden(struct):
                raise ForbiddenStructureException()
            if fails_species_constraints(struct):
                raise ForbiddenStructureException()

        return product_structures

    def is_molecule_forbidden(self, molecule):
        """
        Return ``True`` if the molecule is forbidden in this family, or
        ``False`` otherwise.
        """

        # check family-specific forbidden structures
        if self.forbidden is not None and self.forbidden.is_molecule_forbidden(molecule):
            return True

        # forbid vdw multi-dentate molecules for surface families
        if "surface" in self.label.lower():
            if molecule.get_num_atoms('X') > 1:
                for atom in molecule.atoms:
                    if atom.atomtype.label == 'Xv':
                        return True

        return False

    def _create_reaction(self, reactants, products, is_forward):
        """
        Create and return a new :class:`Reaction` object containing the
        provided `reactants` and `products` as lists of :class:`Molecule`
        objects.
        """

        # Make sure the products are in fact different than the reactants
        if same_species_lists(reactants, products, save_order=self.save_order):
            return None

        # Create and return template reaction object
        reaction = TemplateReaction(
            reactants=reactants if is_forward else products,
            products=products if is_forward else reactants,
            degeneracy=1,
            reversible=self.reversible,
            family=self.label,
            is_forward=is_forward,
            electrons = self.electrons
        )

        if not self.allow_charged_species:
            for spc in (reaction.reactants + reaction.products):
                if spc.get_net_charge() != 0:
                    return None

        if not reaction.is_balanced():
            return None

        # Store the labeled atoms so we can recover them later
        # (e.g. for generating reaction pairs and templates)
        for key, species_list in zip(['reactants', 'products'], [reaction.reactants, reaction.products]):
            for species in species_list:
                reaction.labeled_atoms[key] = dict(reaction.labeled_atoms[key], **species.get_all_labeled_atoms())

        return reaction

    def _match_reactant_to_template(self, reactant, template_reactant):
        """
        Return a complete list of the mappings if the provided reactant
        matches the provided template reactant, or an empty list if not.
        """

        if isinstance(template_reactant, list):
            template_reactant = template_reactant[0]
        if isinstance(template_reactant, Entry):
            struct = template_reactant.item
        else:
            struct = template_reactant

        reactant_contains_surface_site = reactant.contains_surface_site()
        reactant_is_surface_site = reactant.is_surface_site()

        if isinstance(struct, LogicNode):
            mappings = []
            for child_structure in struct.get_possible_structures(self.groups.entries):
                if child_structure.contains_surface_site() != reactant_contains_surface_site:
                    # An adsorbed template can't match a gas-phase species and vice versa
                    continue
                mappings.extend(reactant.find_subgraph_isomorphisms(child_structure, save_order=self.save_order))
            return mappings
        elif isinstance(struct, Group):
            if struct.is_surface_site() != reactant_is_surface_site:
                # An empty surface site group should not match an adsorbate
                return []
            if struct.contains_surface_site() != reactant_contains_surface_site:
                # An adsorbed template can't match a gas-phase species and vice versa
                return []
            return reactant.find_subgraph_isomorphisms(struct, save_order=self.save_order)
        else:
            raise NotImplementedError("Not expecting template of type {}".format(type(struct)))

    def generate_reactions(self, reactants, products=None, prod_resonance=True, delete_labels=True, relabel_atoms=True):
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
                                             Defaults to ``True``, resonance structures are compared.
            delete_labels (bool, optional):  Delete the labeled atoms from each generated reaction (optional).
                                             Default is ``True``, atom labels are deleted.
            relabel_atoms (bool, optional)   Whether to reverse product atom labels of reversible families.
                                             Default is ``True``, atoms are re-labeled.

        Returns:
            List of all reactions containing Molecule objects with the
            specified reactants and products within this family.
            Degenerate reactions are returned as separate reactions.
        """
        reaction_list = []

        # Forward direction (the direction in which kinetics is defined)
        reaction_list.extend(
            self._generate_reactions(reactants=reactants,
                                     products=products,
                                     forward=True,
                                     prod_resonance=prod_resonance,
                                     delete_labels=delete_labels,
                                     relabel_atoms=relabel_atoms,
                                     ))

        if not self.own_reverse and self.reversible:
            # Reverse direction (the direction in which kinetics is not defined)
            reaction_list.extend(
                self._generate_reactions(reactants=reactants,
                                         products=products,
                                         forward=False,
                                         prod_resonance=prod_resonance,
                                         delete_labels=delete_labels,
                                         relabel_atoms=relabel_atoms,
                                         ))
        return reaction_list

    def add_reverse_attribute(self, rxn, react_non_reactive=True):
        """
        For rxn (with species' objects) from families with ownReverse, this method adds a `reverse`
        attribute that contains the reverse reaction information (like degeneracy)

        Returns `True` if successful and `False` if the reverse reaction is forbidden.
        Will raise a `KineticsError` if unsuccessful for other reasons.
        """
        if self.own_reverse and all([spc.has_reactive_molecule() for spc in rxn.products]):
            # Check if the reactants are the same
            same_reactants = 0
            if len(rxn.products) == 2 and rxn.products[0].is_isomorphic(rxn.products[1]):
                same_reactants = 2
            elif len(rxn.products) == 3:
                same_01 = rxn.products[0].is_isomorphic(rxn.products[1])
                same_02 = rxn.products[0].is_isomorphic(rxn.products[2])
                if same_01 and same_02:
                    same_reactants = 3
                elif same_01 or same_02:
                    same_reactants = 2
                elif rxn.products[1].is_isomorphic(rxn.products[2]):
                    same_reactants = 2

            ensure_independent_atom_ids(rxn.products)

            reaction_list = self._generate_reactions([spc.molecule for spc in rxn.products],
                                                     products=rxn.reactants, forward=True,
                                                     react_non_reactive=react_non_reactive)
            reactions = find_degenerate_reactions(reaction_list, same_reactants, kinetics_family=self)
            if len(reactions) == 0:
                logging.error("Expecting one matching reverse reaction, not zero in reaction family {0} for "
                              "forward reaction {1}.\n".format(self.label, str(rxn)))
                logging.error("There is likely a bug in the RMG-database kinetics reaction family involving a "
                              "missing group, missing atomlabels, forbidden groups, etc.")
                for reactant in rxn.reactants:
                    logging.info("Reactant")
                    logging.info(reactant.to_adjacency_list())
                for product in rxn.products:
                    logging.info("Product")
                    logging.info(product.to_adjacency_list())
                logging.error("Debugging why no reaction was found...")
                logging.error("Checking whether the family's forbidden species have affected reaction generation...")
                # Set family's forbidden structures to empty for now to see if reaction gets generated...
                # Note that it is not necessary to check global forbidden structures, because this reaction
                # would not have been formed in the first place.
                temp_object = self.forbidden
                self.forbidden = ForbiddenStructures()  # Initialize with empty one
                try:
                    reaction_list = self._generate_reactions([spc.molecule for spc in rxn.products],
                                                             products=rxn.reactants, forward=True,
                                                             react_non_reactive=react_non_reactive)
                    reactions = find_degenerate_reactions(reaction_list, same_reactants, kinetics_family=self)
                finally:
                    self.forbidden = temp_object
                if (len(reactions) == 1 or
                        (len(reactions) > 1 and
                         all([reactions[0].is_isomorphic(other, check_template_rxn_products=True)
                              for other in reactions]))):
                    logging.error("Error was fixed, the product is a forbidden structure when used as a reactant "
                                  "in the reverse direction.")
                    # This reaction should be forbidden in the forward direction as well
                    return False
                else:
                    logging.error("Still experiencing error: Expecting one matching reverse reaction, not {0} in "
                                  "reaction family {1} for forward reaction {2}."
                                  "\n".format(len(reactions), self.label, str(rxn)))
                    raise KineticsError("Did not find reverse reaction in reaction family {0} for reaction "
                                        "{1}.".format(self.label, str(rxn)))
            elif (len(reactions) > 1 and
                    not all([reactions[0].is_isomorphic(other, strict=False, check_template_rxn_products=True)
                             for other in reactions])):
                logging.error("Expecting one matching reverse reaction. Recieved {0} reactions with "
                              "multiple non-isomorphic ones in reaction family {1} for "
                              "forward reaction {2}.\n".format(len(reactions), self.label, str(rxn)))
                logging.info("Found the following reverse reactions")
                for rxn0 in reactions:
                    logging.info(str(rxn0))
                    for reactant in rxn0.reactants:
                        logging.info("Reactant")
                        logging.info(reactant.to_adjacency_list())
                    for product in rxn0.products:
                        logging.info("Product")
                        logging.info(product.to_adjacency_list())
                raise KineticsError("Found multiple reverse reactions in reaction family {0} for reaction {1}, likely "
                                    "due to inconsistent resonance structure generation".format(self.label, str(rxn)))
            else:
                rxn.reverse = reactions[0]
                return True

    def calculate_degeneracy(self, reaction, resonance=True):
        """
        For a `reaction`  with `Molecule` or `Species` objects given in the direction in which
        the kinetics are defined, compute the reaction-path degeneracy. Can specify whether to consider resonance.

        This method by default adjusts for double counting of identical reactants.
        This should only be adjusted once per reaction. To not adjust for
        identical reactants (since you will be reducing them later in the algorithm), add
        `ignoreSameReactants= True` to this method.
        """
        # Check if the reactants are the same
        # If they refer to the same memory address, then make a deep copy so
        # they can be manipulated independently
        if reaction.is_charge_transfer_reaction():
            # Not implemented yet for charge transfer reactions
            return 1
        reactants = reaction.reactants
        reactants, same_reactants = check_for_same_reactants(reactants)

        # Label reactant atoms for proper degeneracy calculation
        ensure_independent_atom_ids(reactants, resonance=resonance)
        molecule_combos = generate_molecule_combos(reactants)

        reactions = []
        for combo in molecule_combos:
            reactions.extend(self._generate_reactions(combo, products=reaction.products, forward=True,
                                                      prod_resonance=resonance, react_non_reactive=True))

        # remove degenerate reactions
        reactions = find_degenerate_reactions(reactions, same_reactants, template=reaction.template,
                                              kinetics_family=self, resonance=resonance)

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

    def _generate_reactions(self, reactants, products=None, forward=True, prod_resonance=True,
                            react_non_reactive=False, delete_labels=True, relabel_atoms=True):
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
            reactants:          List of Molecules to react.
            products:           List of Molecules or Species of desired product structures (optional).
            forward:            Flag to indicate whether the forward or reverse template should be applied (optional).
                                Default is ``True``, forward template is used.
            prod_resonance:     Flag to generate resonance structures for product checking (optional).
                                Default is ``True``, resonance structures are compared.
            react_non_reactive: Flag to generate reactions between unreactive molecules (optional).
                                Default is ``False``, reactions involving unreactive molecules are not generated.
            delete_labels:      Delete the labeled atoms from each generated reaction (optional).
                                Default is ``True``, atom labels are deleted.
            relabel_atoms (bool, optional)   Whether to reverse product atom labels of reversible families.

        Returns:
            List of all reactions containing Molecule objects with the
                specified reactants and products within this family.
            Degenerate reactions are returned as separate reactions.
        """

        rxn_list = []

        # Wrap each reactant in a list if not already done (this is done to
        # allow for passing multiple resonance structures for each molecule)
        # This also makes a copy of the reactants list so we don't modify the
        # original
        reactants = [reactant if isinstance(reactant, list) else [reactant] for reactant in reactants]

        if forward:
            template = self.forward_template
            reactant_num = self.reactant_num
        elif self.reverse_template is None:
            return []
        else:
            template = self.reverse_template
            reactant_num = self.product_num

        if self.auto_generated and reactant_num != len(reactants):
            return []

        if len(reactants) > len(template.reactants):
            # If the template contains a surface site, we do not want to split it because it will break vdw bonds
            if isinstance(template.reactants[0].item, Group):
                if template.reactants[0].item.contains_surface_site():
                    return []
            # if the family has one template and is bimolecular split template into multiple reactants
            try:
                grps = template.reactants[0].item.split()
                template_reactants = []
                for grp in grps:
                    template_reactants.append(grp)
            except AttributeError:
                template_reactants = [x.item for x in template.reactants]
        else:
            template_reactants = [x.item for x in template.reactants]

        # Unimolecular reactants: A --> products
        if len(reactants) == 1 and len(template_reactants) == 1:

            # Iterate over all resonance isomers of the reactant
            for molecule in reactants[0]:
                if molecule.reactive or react_non_reactive:  # don't react non representative resonance isomers unless
                    # explicitly desired (e.g., when called from calculate_degeneracy)
                    mappings = self._match_reactant_to_template(molecule, template_reactants[0])
                    for mapping in mappings:
                        reactant_structures = [molecule]
                        try:
                            product_structures = self._generate_product_structures(reactant_structures,
                                                                                   [mapping],
                                                                                   forward,
                                                                                   relabel_atoms)
                        except ForbiddenStructureException:
                            pass
                        else:
                            if product_structures is not None:
                                rxn = self._create_reaction(reactant_structures, product_structures, forward)
                                if rxn:
                                    rxn_list.append(rxn)
        # Bimolecular reactants: A + B --> products
        elif len(reactants) == 2 and len(template_reactants) == 2:

            molecules_a = reactants[0]
            molecules_b = reactants[1]

            if 'adsorption' in self.label.lower() and forward:
                if molecules_a[0].contains_surface_site() and molecules_b[0].contains_surface_site():
                    if 'vdw' in self.label.lower():
                        # can adsorb vdW species to the surface, so continue on
                        pass
                    else:
                        # Can't adsorb something that's already adsorbed.
                        # Both reactants either contain or are a surface site.
                        return []

            # Iterate over all resonance isomers of the reactant
            for molecule_a in molecules_a:
                for molecule_b in molecules_b:
                    if (molecule_a.reactive and molecule_b.reactive) or react_non_reactive:

                        # Reactants stored as A + B
                        mappings_a = self._match_reactant_to_template(molecule_a, template_reactants[0])
                        mappings_b = self._match_reactant_to_template(molecule_b, template_reactants[1])

                        # Iterate over each pair of matches (A, B)
                        for map_a in mappings_a:
                            for map_b in mappings_b:
                                # Reverse the order of reactants in case we have a family with only one reactant tree
                                # that can produce different products depending on the order of reactants
                                reactant_structures = [molecule_b, molecule_a]
                                try:
                                    product_structures = self._generate_product_structures(reactant_structures,
                                                                                           [map_b, map_a],
                                                                                           forward,
                                                                                           relabel_atoms)
                                except ForbiddenStructureException:
                                    pass
                                else:
                                    if product_structures is not None:
                                        rxn = self._create_reaction(reactant_structures, product_structures, forward)
                                        if rxn:
                                            rxn_list.append(rxn)

                        # Only check for swapped reactants if they are different
                        if reactants[0] is not reactants[1]:

                            # Reactants stored as B + A
                            mappings_a = self._match_reactant_to_template(molecule_a, template_reactants[1])
                            mappings_b = self._match_reactant_to_template(molecule_b, template_reactants[0])

                            # Iterate over each pair of matches (A, B)
                            for map_a in mappings_a:
                                for map_b in mappings_b:
                                    reactant_structures = [molecule_a, molecule_b]
                                    try:
                                        product_structures = self._generate_product_structures(reactant_structures,
                                                                                               [map_a, map_b],
                                                                                               forward,
                                                                                               relabel_atoms)
                                    except ForbiddenStructureException:
                                        pass
                                    else:
                                        if product_structures is not None:
                                            rxn = self._create_reaction(reactant_structures, product_structures,
                                                                        forward)
                                            if rxn:
                                                rxn_list.append(rxn)

        # Termolecular reactants: A + B + C --> products
        elif len(reactants) == 2 and len(template_reactants) == 3:
            """
            Two reactants but a termolecular template.
            Could be A + X + X <=> BX + CX (dissociative adsorption)
            or A + X + X <=> AXX (bidentate adsorption)
            in which case, if one of the two reactants is an X
            then we have a match and can just use it twice.
            """
            template_sites = [r for r in template_reactants if r.is_surface_site()]
            if len(template_sites) == 2:
                # Two surface sites in template. If there's a site in the reactants, use it twice.
                if reactants[0][0].is_surface_site() and not reactants[1][0].is_surface_site():
                    site1 = reactants[0][0]
                    site2 = deepcopy(reactants[0][0])
                    adsorbate_molecules = reactants[1]
                    reactants.append([site2])
                elif reactants[1][0].is_surface_site() and not reactants[0][0].is_surface_site():
                    site1 = reactants[1][0]
                    site2 = deepcopy(reactants[1][0])
                    adsorbate_molecules = reactants[0]
                    reactants.append([site2])
                else:
                    # No reaction with these reactants in this template
                    return []

                for r in template_reactants:
                    if not r.is_surface_site():
                        template_adsorbate = r
                        break
                else:
                    raise KineticsError("Couldn't find non-site in template {0!r}".format(template))

                mappings_a = self._match_reactant_to_template(site1, template_sites[0])
                mappings_b = self._match_reactant_to_template(site2, template_sites[1])
                for adsorbateMolecule in adsorbate_molecules:
                    mappings_c = self._match_reactant_to_template(adsorbateMolecule, template_adsorbate)
                    for map_a, map_b, map_c in itertools.product(mappings_a, mappings_b, mappings_c):
                        reactant_structures = [site1, site2, adsorbateMolecule]
                        # should be in same order as reaction template recipe?
                        try:
                            product_structures = self._generate_product_structures(reactant_structures,
                                                                                   [map_a, map_b, map_c],
                                                                                   forward,
                                                                                   relabel_atoms)
                        except ForbiddenStructureException:
                            pass
                        else:
                            if product_structures is not None:
                                rxn = self._create_reaction(reactant_structures, product_structures, forward)
                                if rxn:
                                    rxn_list.append(rxn)
            else:
                # _generate_reactions was called with mismatched number of reactants and templates
                return []

        elif len(reactants) == 3 and len(template_reactants) == 3:
            """
            This could be a surface reaction
                A + X + X <=> BX + CX    (dissociative adsorption)
                A + X + X <=> AXX        (bidentate adsorption)
                ABX + X + X <=> AXX + BX (dissociation to bidentate)
            or a termolecular gas phase reaction
                A + B + C <=> stuff
            We check the two scenarios in that order.
            """
            template_sites = [r for r in template_reactants if r.is_surface_site()]
            if len(template_sites) == 2:
                """
                Three reactants and a termolecular template.
                Could be A + X + X <=> BX + CX (dissociative adsorption)
                or A + X + X <=> AXX (bidentate adsorption)
                that was first found in the reverse direction
                and so is being passed in with all three reactants identified.
                """
                # Should be 2 surface sites in reactants too.
                # Find them, and find mappings of the other
                m1, m2, m3 = (r[0] for r in reactants)
                if m1.is_surface_site() and m2.is_surface_site() and not m3.is_surface_site():
                    site1, site2 = m1, m2
                    adsorbate_molecules = reactants[2]
                elif m1.is_surface_site() and not m2.is_surface_site() and m3.is_surface_site():
                    site1, site2 = m1, m3
                    adsorbate_molecules = reactants[1]
                elif not m1.is_surface_site() and m2.is_surface_site() and m3.is_surface_site():
                    site1, site2 = m2, m3
                    adsorbate_molecules = reactants[0]
                else:
                    # Three reactants not containing two surface sites
                    return []

                for r in template_reactants:
                    if not r.is_surface_site():
                        template_adsorbate = r
                        break
                else:
                    raise KineticsError("Couldn't find non-site in template {0!r}".format(template))

                mappings_a = self._match_reactant_to_template(site1, template_sites[0])
                mappings_b = self._match_reactant_to_template(site2, template_sites[1])
                for adsorbateMolecule in adsorbate_molecules:
                    mappings_c = self._match_reactant_to_template(adsorbateMolecule, template_adsorbate)
                    # this just copied/pasted from above - not checked
                    for map_a, map_b, map_c in itertools.product(mappings_a, mappings_b, mappings_c):
                        reactant_structures = [site1, site2, adsorbateMolecule]
                        try:
                            product_structures = self._generate_product_structures(reactant_structures,
                                                                                   [map_a, map_b, map_c],
                                                                                   forward,
                                                                                   relabel_atoms)
                        except ForbiddenStructureException:
                            pass
                        else:
                            if product_structures is not None:
                                rxn = self._create_reaction(reactant_structures, product_structures, forward)
                                if rxn:
                                    rxn_list.append(rxn)

            else:
                """
                Not a bidentate surface reaction, just a gas-phase
                Trimolecular reactants: A + B + C --> products
                """
                molecules_a = reactants[0]
                molecules_b = reactants[1]
                molecules_c = reactants[2]

                # Iterate over all resonance isomers of the reactants
                for molecule_a in molecules_a:
                    for molecule_b in molecules_b:
                        for molecule_c in molecules_c:

                            def generate_products_and_reactions(order):
                                """
                                order = (0, 1, 2) corresponds to reactants stored as A + B + C, etc.
                                """
                                _mappings_a = self._match_reactant_to_template(molecule_a, template_reactants[order[0]])
                                _mappings_b = self._match_reactant_to_template(molecule_b, template_reactants[order[1]])
                                _mappings_c = self._match_reactant_to_template(molecule_c, template_reactants[order[2]])

                                # Iterate over each pair of matches (A, B, C)
                                for _map_a in _mappings_a:
                                    for _map_b in _mappings_b:
                                        for _map_c in _mappings_c:
                                            _reactantStructures = [molecule_a, molecule_b, molecule_c]
                                            _maps = [_map_a, _map_b, _map_c]
                                            # Reorder reactants in case we have a family with fewer reactant trees than
                                            # reactants and different reactant orders can produce different products
                                            _reactantStructures = [_reactantStructures[_i] for _i in order]
                                            _maps = [_maps[_i] for _i in order]
                                            try:
                                                _productStructures = self._generate_product_structures(
                                                    _reactantStructures,
                                                    _maps,
                                                    forward,
                                                    relabel_atoms)
                                            except ForbiddenStructureException:
                                                pass
                                            else:
                                                if _productStructures is not None:
                                                    _rxn = self._create_reaction(_reactantStructures,
                                                                                 _productStructures,
                                                                                 forward)
                                                    if _rxn:
                                                        rxn_list.append(_rxn)

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

        # ToDo: try to remove this hard-coding of reaction family name..
        if not forward and ('adsorption' in self.label.lower() or 'eleyrideal' in self.label.lower()):
            # Desorption should have desorbed something (else it was probably bidentate)
            # so delete reactions that don't make a gas-phase desorbed product
            # Eley-Rideal reactions should have one gas-phase product in the reverse direction

            # Determine how many surf reactants we expect based on the template
            n_surf_expected = len([r for r in self.forward_template.reactants if r.item.contains_surface_site()])

            # Now iterate through the reactions and toss them out if the number of surface reactants
            # does not match the expected number
            pruned_list = []
            for reaction in rxn_list:
                n_surf_reaction = len([r for r in reaction.reactants if r.contains_surface_site()])
                if n_surf_expected != n_surf_reaction:
                    logging.debug("Removing {0} reaction {1!s} with no desorbed species".format(self.label, reaction))
                    continue
                else:
                    pruned_list.append(reaction)
            rxn_list = pruned_list

        # If products is given, remove reactions from the reaction list that
        # don't generate the given products
        if products is not None:
            rxn_list0 = rxn_list[:]
            rxn_list = []
            for reaction in rxn_list0:
                products0 = reaction.products if forward else reaction.reactants
                # Only keep reactions which give the requested products
                # If prod_resonance=True, then use strict=False to consider all resonance structures
                if same_species_lists(products, products0, strict=not prod_resonance, save_order=self.save_order):
                    rxn_list.append(reaction)

        # Determine the reactant-product pairs to use for flux analysis
        # Also store the reaction template (useful so we can easily get the kinetics later)
        for reaction in rxn_list:

            # Restore the labeled atoms long enough to generate some metadata
            for reactant in reaction.reactants:
                reactant.clear_labeled_atoms()
            for label, atom in reaction.labeled_atoms['reactants'].items():
                if isinstance(atom, list):
                    for atm in atom:
                        atm.label = label
                else:
                    atom.label = label

            # Generate metadata about the reaction that we will need later
            reaction.pairs = self.get_reaction_pairs(reaction)
            reaction.template = self.get_reaction_template_labels(reaction)

            if delete_labels:
                # Unlabel the atoms for both reactants and products
                for species in itertools.chain(reaction.reactants, reaction.products):
                    species.clear_labeled_atoms()

                # We're done with the labeled atoms, so delete the attribute
                del reaction.labeled_atoms

            # Mark reaction reversibility
            reaction.reversible = self.reversible

        # This reaction list has only checked for duplicates within itself, not
        # with the global list of reactions
        return rxn_list

    def get_reaction_pairs(self, reaction):
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
                    pairs.append([reactant, product])
        elif self.label.lower() in ('h_abstraction','f_abstraction','cl_abstraction','br_abstraction'):
            # Hardcoding for hydrogen abstraction: pair the reactant containing
            # *1 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].contains_labeled_atom('*1'):
                if reaction.products[0].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[1].contains_labeled_atom('*1'):
                if reaction.products[1].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[0].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
        elif self.label.lower() in ['disproportionation', 'co_disproportionation', 'korcek_step1_cat']:
            # Hardcoding for disproportionation, co_disproportionation, korcek_step1_cat:
            # pair the reactant containing *1 with the product containing *1
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].contains_labeled_atom('*1'):
                if reaction.products[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[1].contains_labeled_atom('*1'):
                if reaction.products[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
        elif self.label.lower() in ['substitution_o', 'substitutions']:
            # Hardcoding for Substitution_O: pair the reactant containing
            # *2 with the product containing *3 and vice versa
            assert len(reaction.reactants) == len(reaction.products) == 2
            if reaction.reactants[0].contains_labeled_atom('*2'):
                if reaction.products[0].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[1].contains_labeled_atom('*2'):
                if reaction.products[1].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[0].contains_labeled_atom('*3'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
        elif self.label.lower() == 'baeyer-villiger_step1_cat':
            # Hardcoding for Baeyer-Villiger_step1_cat: pair the two reactants
            # with the Criegee intermediate and pair the catalyst with itself
            assert len(reaction.reactants) == 3 and len(reaction.products) == 2
            if reaction.reactants[0].contains_labeled_atom('*5'):
                if reaction.products[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                elif reaction.products[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[0]])
            elif reaction.reactants[1].contains_labeled_atom('*5'):
                if reaction.products[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.products[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
            elif reaction.reactants[2].contains_labeled_atom('*5'):
                if reaction.products[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[2], reaction.products[1]])
                elif reaction.products[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[2], reaction.products[0]])
        elif self.label.lower() == 'baeyer-villiger_step2_cat':
            # Hardcoding for Baeyer-Villiger_step2_cat: pair the Criegee
            # intermediate with the two products and the catalyst with itself
            assert len(reaction.reactants) == 2 and len(reaction.products) == 3
            if reaction.products[0].contains_labeled_atom('*7'):
                if reaction.reactants[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                elif reaction.reactants[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                    pairs.append([reaction.reactants[0], reaction.products[0]])
            elif reaction.products[1].contains_labeled_atom('*7'):
                if reaction.reactants[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                elif reaction.reactants[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
            elif reaction.products[2].contains_labeled_atom('*7'):
                if reaction.reactants[0].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[0], reaction.products[0]])
                    pairs.append([reaction.reactants[0], reaction.products[1]])
                    pairs.append([reaction.reactants[1], reaction.products[2]])
                elif reaction.reactants[1].contains_labeled_atom('*1'):
                    pairs.append([reaction.reactants[1], reaction.products[0]])
                    pairs.append([reaction.reactants[1], reaction.products[1]])
                    pairs.append([reaction.reactants[0], reaction.products[2]])
        elif reaction.is_surface_reaction():
            # remove vacant active sites from consideration
            reactants = [sp for sp in reaction.reactants if not sp.is_surface_site()]
            products = [sp for sp in reaction.products if not sp.is_surface_site()]
            if len(reactants) == 1 or len(products) == 1:
                # When there is only one reactant (or one product), it is paired
                # with each of the products (reactants)
                for reactant in reactants:
                    for product in products:
                        pairs.append([reactant, product])
            elif self.label.lower() == 'surface_abstraction':
                # Hardcoding for surface abstraction: pair the reactant containing
                # *1 with the product containing *3 and vice versa
                assert len(reaction.reactants) == len(reaction.products) == 2
                if reaction.reactants[0].contains_labeled_atom('*1'):
                    if reaction.products[0].contains_labeled_atom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[0]])
                        pairs.append([reaction.reactants[1], reaction.products[1]])
                    elif reaction.products[1].contains_labeled_atom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[1]])
                        pairs.append([reaction.reactants[1], reaction.products[0]])
                elif reaction.reactants[1].contains_labeled_atom('*1'):
                    if reaction.products[1].contains_labeled_atom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[0]])
                        pairs.append([reaction.reactants[1], reaction.products[1]])
                    elif reaction.products[0].contains_labeled_atom('*3'):
                        pairs.append([reaction.reactants[0], reaction.products[1]])
                        pairs.append([reaction.reactants[1], reaction.products[0]])
        if not pairs:
            logging.debug('Preset mapping missing for determining reaction pairs for family {0!s}, '
                          'falling back to Reaction.generate_pairs'.format(self.label))

        return pairs

    def get_reaction_template(self, reaction):
        """
        For a given `reaction` with properly-labeled :class:`Molecule` objects
        as the reactants, determine the most specific nodes in the tree that
        describe the reaction.
        """
        return self.groups.get_reaction_template(reaction)

    def get_kinetics_for_template(self, template, degeneracy=1, method='rate rules'):
        """
        Return an estimate of the kinetics for a reaction with the given
        `template` and reaction-path `degeneracy`. There is currently only one method to use:
        'rate rules' (old RMG-Java behavior, and default RMG-Py behavior). Group additivity was removed in August 2023.
        
        Returns a tuple (kinetics, entry):
        If it's estimated via 'rate rules' and an exact match is found in the tree,
        then the entry is returned as the second element of the tuple.
        But if an average is used, then the tuple returned is (kinetics, None).

        """
        if method.lower() == 'rate rules':
            return self.estimate_kinetics_using_rate_rules(template, degeneracy)  # This returns kinetics and entry data
        else:
            raise ValueError('Invalid value "{0}" for method parameter; '
                             'currently only "rate rules" is supported.'.format(method))

    def get_kinetics_from_depository(self, depository, reaction, template, degeneracy):
        """
        Search the given `depository` in this kinetics family for kinetics
        for the given `reaction`. Returns a list of all of the matching
        kinetics, the corresponding entries, and ``True`` if the kinetics
        match the forward direction or ``False`` if they match the reverse
        direction.
        """
        kinetics_list = []
        entries = depository.entries.values()
        for entry in entries:
            if entry.item.is_isomorphic(reaction):
                kinetics_list.append(
                    [deepcopy(entry.data), entry, entry.item.is_isomorphic(reaction, either_direction=False)])
        for kinetics, entry, is_forward in kinetics_list:
            if kinetics is not None:
                kinetics.comment += "Matched reaction {0} {1} in {2}\nThis reaction matched rate rule {3}".format(
                    entry.index,
                    entry.label,
                    depository.label,
                    '[{0}]'.format(';'.join([g.label for g in template])))
                kinetics.comment += "\nfamily: {}".format(self.label)
                if entry.metal:
                    kinetics.comment += "\nmetal: {}".format(self.metal)
                if entry.facet:
                    kinetics.comment += "\nfacet: {}".format(self.facet)
                if entry.facet:
                    kinetics.comment += "\nsite: {}".format(self.site)
        return kinetics_list

    def _select_best_kinetics(self, kinetics_list):
        """
        For a given set of kinetics `kinetics_list`, return the kinetics deemed
        to be the "best". This is determined to be the one with the lowest
        non-zero rank that occurs first (has the lowest index).
        """
        if any([x[1].rank == 0 for x in kinetics_list]) and not all([x[1].rank == 0 for x in kinetics_list]):
            kinetics_list = [x for x in kinetics_list if x[1].rank != 0]
        kinetics_list.sort(key=lambda x: (x[1].rank, x[1].index))
        return kinetics_list[0]

    def get_kinetics(self, reaction, template_labels, degeneracy=1, estimator='', return_all_kinetics=True):
        """
        Return the kinetics for the given `reaction` by searching the various
        depositories as well as generating a result using the user-specified `estimator`.
        Currently, only 'rate rules' is a supported estimator.  Unlike
        the regular :meth:`get_kinetics()` method, this returns a list of
        results, with each result comprising of

        1. the kinetics
        2. the source - this will be `None` if from a template estimate
        3. the entry  - this will be `None` if from a template estimate
        4. is_forward a boolean denoting whether the matched entry is in the same
           direction as the inputted reaction. This will always be True if using
           rates rules. This can be `True` or `False` if using
           a depository

        If return_all_kinetics==False, only the first (best?) matching kinetics is returned.
        """
        kinetics_list = []

        depositories = self.depositories[:]

        template = self.retrieve_template(template_labels)

        # Check the various depositories for kinetics
        for depository in depositories:
            kinetics_list0 = self.get_kinetics_from_depository(depository, reaction, template, degeneracy)
            if len(kinetics_list0) > 0 and not return_all_kinetics:
                kinetics, entry, is_forward = self._select_best_kinetics(kinetics_list0)
                return kinetics, depository, entry, is_forward
            else:
                for kinetics, entry, is_forward in kinetics_list0:
                    kinetics_list.append([kinetics, depository, entry, is_forward])

        # If estimator type of rate rules is given, retrieve the kinetics. 
        # TODO: Since group additivity was removed, this logic can be condensed into just 1 branch.

        if estimator:
            try:
                kinetics, entry = self.get_kinetics_for_template(template, degeneracy, method=estimator)
            except Exception:
                logging.error("Error getting kinetics for reaction {0!s}.\n{0!r}".format(reaction))
                raise

            if kinetics:
                if not return_all_kinetics:
                    return kinetics, estimator, entry, True
                kinetics_list.append([kinetics, estimator, entry, True])
        # If no estimation method was given, prioritize rate rule estimation. 
        else:
            try:
                kinetics, entry = self.get_kinetics_for_template(template, degeneracy, method='rate rules')
                if not return_all_kinetics:
                    return kinetics, 'rate rules', entry, True
                kinetics_list.append([kinetics, 'rate rules', entry, True])
            except KineticsError:
                # If kinetics were undeterminable for rate rules estimation, do nothing.
                pass

        if not return_all_kinetics:
            raise UndeterminableKineticsError(reaction)

        return kinetics_list

    def estimate_kinetics_using_rate_rules(self, template, degeneracy=1):
        """
        Determine the appropriate kinetics for a reaction with the given
        `template` using rate rules.

        Returns a tuple (kinetics, entry) where `entry` is the database
        entry used to determine the kinetics only if it is an exact match,
        and is None if some averaging or use of a parent node took place.
        """
        kinetics, entry = self.rules.estimate_kinetics(template, degeneracy)

        return kinetics, entry

    def get_reaction_template_labels(self, reaction):
        """
        Retrieve the template for the reaction and
        return the corresponding labels for each of the
        groups in the template.
        """
        template = self.get_reaction_template(reaction)

        template_labels = []
        for entry in template:
            template_labels.append(entry.label)

        return template_labels

    def retrieve_template(self, template_labels):
        """
        Reconstruct the groups associated with the
        labels of the reaction template and
        return a list.
        """
        template = []
        for label in template_labels:
            template.append(self.groups.entries[label])

        return template

    def get_labeled_reactants_and_products(self, reactants, products, relabel_atoms=True):
        """
        Given `reactants`, a list of :class:`Molecule` objects, and products, a list of
        :class:`Molecule` objects, return two new lists of :class:`Molecule` objects with
        atoms labeled: one for reactants, one for products. Returned molecules are totally
        new entities in memory so input molecules `reactants` and `products` won't be affected.
        If RMG cannot find appropriate labels, (None, None) will be returned.
        If ``relabel_atoms`` is ``True``, product atom labels of reversible families
        will be reversed to assist in identifying forbidden structures.
        """
        template = self.forward_template
        reactants0 = [reactant.copy(deep=True) for reactant in reactants]

        if self.auto_generated and self.reactant_num != len(reactants):
            return None, None

        if len(reactants) > len(template.reactants):
            # if the family has one template and is bimolecular split template into multiple reactants
            try:
                grps = template.reactants[0].item.split()
                template_reactants = []
                for grp in grps:
                    template_reactants.append(grp)
            except AttributeError:
                template_reactants = [x.item for x in template.reactants]
        else:
            template_reactants = [x.item for x in template.reactants]

        if len(reactants0) == 1:
            molecule = reactants0[0]
            mappings = self._match_reactant_to_template(molecule, template_reactants[0])
            mappings = [[map0] for map0 in mappings]
            num_mappings = len(mappings)
            reactant_structures = [molecule]
        elif len(reactants0) == 2:
            molecule_a = reactants0[0]
            molecule_b = reactants0[1]
            # get mappings in forward direction
            mappings_a = self._match_reactant_to_template(molecule_a, template_reactants[0])
            mappings_b = self._match_reactant_to_template(molecule_b, template_reactants[1])
            mappings = list(itertools.product(mappings_a, mappings_b))
            # get mappings in the reverse direction
            mappings_a = self._match_reactant_to_template(molecule_a, template_reactants[1])
            mappings_b = self._match_reactant_to_template(molecule_b, template_reactants[0])
            mappings.extend(list(itertools.product(mappings_a, mappings_b)))

            reactant_structures = [molecule_a, molecule_b]
            num_mappings = len(mappings_a) * len(mappings_b)
        elif len(reactants0) == 3:
            molecule_a = reactants0[0]
            molecule_b = reactants0[1]
            molecule_c = reactants0[2]
            # Get mappings for all permutations of reactants
            mappings = []
            for order in itertools.permutations(range(3), 3):
                mappings_a = self._match_reactant_to_template(molecule_a, template_reactants[order[0]])
                mappings_b = self._match_reactant_to_template(molecule_b, template_reactants[order[1]])
                mappings_c = self._match_reactant_to_template(molecule_c, template_reactants[order[2]])
                mappings.extend(list(itertools.product(mappings_a, mappings_b, mappings_c)))

            reactant_structures = [molecule_a, molecule_b, molecule_c]
            num_mappings = len(mappings_a) * len(mappings_b) * len(mappings_c)
        else:
            raise IndexError('You have {0} reactants, which is unexpected!'.format(len(reactants)))

        for mapping in mappings:
            try:
                product_structures = self._generate_product_structures(reactant_structures, mapping, forward=True, relabel_atoms=relabel_atoms)
            except ForbiddenStructureException:
                pass
            else:
                if product_structures is not None:
                    if same_species_lists(list(products), list(product_structures), save_order=self.save_order):
                        return reactant_structures, product_structures
                    else:
                        continue

        # if there're some mapping available but cannot match the provided products
        # raise exception
        if num_mappings > 0:
            raise ActionError('Something wrong with products that RMG cannot find a match!')

        return None, None

    def add_atom_labels_for_reaction(self, reaction, output_with_resonance=True, save_order=False, relabel_atoms=False):
        """
        Apply atom labels on a reaction using the appropriate atom labels from
        this reaction family.

        The reaction is modified in place containing species objects with the
        atoms labeled. If output_with_resonance is True, all resonance structures
        are generated with labels. If false, only the first resonance structure
        successfully able to map to the reaction is used. None is returned.
        If ``save_order`` is ``True`` the atom order is reset after performing atom isomorphism.
        If ``relabel_atoms`` is ``True``, product atom labels of reversible families
        will be reversed to assist in identifying forbidden structures.
        """
        # make sure we start with reaction with species objects
        reaction.ensure_species(reactant_resonance=False, product_resonance=False, save_order=save_order)

        reactants = reaction.reactants
        products = reaction.products
        # ensure all species are independent references
        if len(reactants + products) > len(set([id(s) for s in reactants + products])):
            logging.debug('Copying reactants and products for reaction {} since they have '
                          'identical species references'.format(reaction))
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
                labeled_reactants, labeled_products = self.get_labeled_reactants_and_products(reactant_pair, product_pair, relabel_atoms=relabel_atoms)
                if labeled_reactants is not None:
                    break
            except ActionError:
                # must have gotten the wrong pair
                pass
        if labeled_reactants is None or labeled_products is None:
            raise ActionError("Could not find labeled reactants for reaction {} from "
                              "family {}.".format(reaction, self.label))

        # place the molecules in reaction's species object
        # this prevents overwriting of attributes of species objects by this method
        for index, species in enumerate(products):
            for labeled_molecule in labeled_products:
                if species.is_isomorphic(labeled_molecule, save_order=self.save_order):
                    species.molecule = [labeled_molecule]
                    reaction.products[index] = species
                    labeled_products.remove(labeled_molecule)
                    break
            else:
                raise ActionError('Could not find isomorphic molecule to fit the original product {} from '
                                  'reaction {}'.format(species, reaction))
        for index, species in enumerate(reactants):
            for labeled_molecule in labeled_reactants:
                if species.is_isomorphic(labeled_molecule, save_order=self.save_order):
                    species.molecule = [labeled_molecule]
                    reaction.reactants[index] = species
                    labeled_reactants.remove(labeled_molecule)
                    break
            else:
                raise ActionError('Could not find isomorphic molecule to fit the original reactant {} from '
                                  'reaction {}'.format(species, reaction))

        if output_with_resonance:
            # convert the molecules to species objects with resonance structures
            for species in reaction.reactants + reaction.products:
                species.generate_resonance_structures()

    def get_training_depository(self):
        """
        Returns the `training` depository from self.depositories
        """
        for depository in self.depositories:
            if depository.label.endswith('training'):
                return depository
        else:
            raise DatabaseError('Could not find training depository in family {0}.'.format(self.label))

    def add_entry(self, parent, grp, name):
        """
        Adds a group entry with parent parent
        group structure grp
        and group name name
        """
        ind = len(self.groups.entries) - 1
        entry = Entry(index=ind, label=name, item=grp, parent=parent)
        self.groups.entries[name] = entry
        self.rules.entries[name] = []
        if entry.parent:
            entry.parent.children.append(entry)

    def _split_reactions(self, rxns, newgrp):
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
        new_inds = []

        for i, rxn in enumerate(rxns):
            rmol = rxn.reactants[0].molecule[0]
            for r in rxn.reactants[1:]:
                rmol = rmol.merge(r.molecule[0])

            rmol.identify_ring_membership()

            if rmol.is_subgraph_isomorphic(newgrp, generate_initial_map=True, save_order=True):
                new.append(rxn)
                new_inds.append(i)
            else:
                comp.append(rxn)

        return new, comp, new_inds

    def reaction_matches(self, rxn, grp):
        rmol = rxn.reactants[0].molecule[0]
        for r in rxn.reactants[1:]:
            rmol = rmol.merge(r.molecule[0])
        rmol.identify_ring_membership()
        return rmol.is_subgraph_isomorphic(grp, generate_initial_map=True, save_order=True)

    def eval_ext(self, parent, ext, extname, template_rxn_map, obj=None, T=1000.0):
        """
        evaluates the objective function obj
        for the extension ext with name extname to the parent entry parent
        """
        rxns = template_rxn_map[parent.label]
        new, old, new_inds = self._split_reactions(rxns, ext)
        if len(new) == 0:
            return np.inf, False
        elif len(old) == 0:
            return np.inf, True
        else:
            if obj:
                ob, boo = get_objective_function(new, old, obj, T=T)
            else:
                ob, boo = get_objective_function(new, old, T=T)
            return ob, True

    def get_extension_edge(self, parent, template_rxn_map, obj, T, iter_max=np.inf, iter_item_cap=np.inf):
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
        out_exts = [[]]
        grps = [[parent.item]]
        names = [parent.label]
        first_time = True
        gave_up_split = False

        n_splits = len(template_rxn_map[parent.label][0].reactants)
        iter = 0

        while grps[iter] != []:
            grp = grps[iter][-1]

            exts = grp.get_extensions(basename=names[-1], n_splits=n_splits)

            reg_dict = dict()
            ext_inds = []
            for i, (grp2, grpc, name, typ, indc) in enumerate(exts):

                if typ != 'intNewBondExt' and typ != 'extNewBondExt' and (typ, indc) not in reg_dict.keys():
                    # first list is all extensions that match at least one reaction
                    # second is extensions that match all reactions
                    reg_dict[(typ, indc)] = ([], [])
                val, boo = self.eval_ext(parent, grp2, name, template_rxn_map, obj, T)

                if val != np.inf:
                    out_exts[-1].append(exts[i])  # this extension splits reactions (optimization dim)
                    if typ == 'atomExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].atomtype)
                    elif typ == 'elExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].radical_electrons)
                    elif typ == 'bondExt':
                        reg_dict[(typ, indc)][0].extend(grp2.get_bond(grp2.atoms[indc[0]], grp2.atoms[indc[1]]).order)

                elif boo:  # this extension matches all reactions (regularization dim)
                    if typ == 'intNewBondExt' or typ == 'extNewBondExt':
                        # these are bond formation extensions, we want to expand these until we get splits
                        ext_inds.append(i)
                    elif typ == 'atomExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].atomtype)
                        reg_dict[(typ, indc)][1].extend(grp2.atoms[indc[0]].atomtype)
                    elif typ == 'elExt':
                        reg_dict[(typ, indc)][0].extend(grp2.atoms[indc[0]].radical_electrons)
                        reg_dict[(typ, indc)][1].extend(grp2.atoms[indc[0]].radical_electrons)
                    elif typ == 'bondExt':
                        reg_dict[(typ, indc)][0].extend(grp2.get_bond(grp2.atoms[indc[0]], grp2.atoms[indc[1]]).order)
                        reg_dict[(typ, indc)][1].extend(grp2.get_bond(grp2.atoms[indc[0]], grp2.atoms[indc[1]]).order)
                    elif typ == 'ringExt':
                        reg_dict[(typ, indc)][1].append(True)
                else:
                    # this extension matches no reactions
                    if typ == 'ringExt':
                        reg_dict[(typ, indc)][0].append(False)
                        reg_dict[(typ, indc)][1].append(False)

            for typr, indcr in reg_dict.keys():  # have to label the regularization dimensions in all relevant groups
                reg_val = reg_dict[(typr, indcr)]

                if first_time and parent.children == []:
                    # parent
                    if typr != 'intNewBondExt' and typr != 'extNewBondExt':  # these dimensions should be regularized
                        if typr == 'atomExt':
                            grp.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                        elif typr == 'elExt':
                            grp.atoms[indcr[0]].reg_dim_u = list(reg_val)
                        elif typr == 'ringExt':
                            grp.atoms[indcr[0]].reg_dim_r = list(reg_val)
                        elif typr == 'bondExt':
                            atms = grp.atoms
                            bd = grp.get_bond(atms[indcr[0]], atms[indcr[1]])
                            bd.reg_dim = list(reg_val)

                # extensions being sent out
                if typr != 'intNewBondExt' and typr != 'extNewBondExt':  # these dimensions should be regularized
                    for grp2, grpc, name, typ, indc in out_exts[-1]:  # returned groups
                        if typr == 'atomExt':
                            grp2.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                        elif typr == 'elExt':
                            grp2.atoms[indcr[0]].reg_dim_u = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_u = list(reg_val)
                        elif typr == 'ringExt':
                            grp2.atoms[indcr[0]].reg_dim_r = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_r = list(reg_val)
                        elif typr == 'bondExt':
                            atms = grp2.atoms
                            bd = grp2.get_bond(atms[indcr[0]], atms[indcr[1]])
                            bd.reg_dim = [list(set(bd.order) & set(reg_val[0])), list(set(bd.order) & set(reg_val[1]))]
                            if grpc:
                                atms = grpc.atoms
                                bd = grpc.get_bond(atms[indcr[0]], atms[indcr[1]])
                                bd.reg_dim = [list(set(bd.order) & set(reg_val[0])),
                                              list(set(bd.order) & set(reg_val[1]))]

            # extensions being expanded
            for typr, indcr in reg_dict.keys():  # have to label the regularization dimensions in all relevant groups
                reg_val = reg_dict[(typr, indcr)]
                if typr != 'intNewBondExt' and typr != 'extNewBondExt':  # these dimensions should be regularized
                    for ind2 in ext_inds:  # groups for expansion
                        grp2, grpc, name, typ, indc = exts[ind2]
                        if typr == 'atomExt':
                            grp2.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_atm = list(reg_val)
                        elif typr == 'elExt':
                            grp2.atoms[indcr[0]].reg_dim_u = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_u = list(reg_val)
                        elif typr == 'ringExt':
                            grp2.atoms[indcr[0]].reg_dim_r = list(reg_val)
                            if grpc:
                                grpc.atoms[indcr[0]].reg_dim_r = list(reg_val)
                        elif typr == 'bondExt':
                            atms = grp2.atoms
                            bd = grp2.get_bond(atms[indcr[0]], atms[indcr[1]])
                            bd.reg_dim = [list(set(bd.order) & set(reg_val[0])), list(set(bd.order) & set(reg_val[1]))]
                            if grpc:
                                atms = grpc.atoms
                                bd = grpc.get_bond(atms[indcr[0]], atms[indcr[1]])
                                bd.reg_dim = [list(set(bd.order) & set(reg_val[0])),
                                              list(set(bd.order) & set(reg_val[1]))]

            out_exts.append([])
            grps[iter].pop()
            names.pop()

            for ind in ext_inds:  # collect the groups to be expanded
                grpr, grpcr, namer, typr, indcr = exts[ind]
                if len(grps) == iter+1:
                    grps.append([])
                grps[iter+1].append(grpr)
                names.append(namer)

            if first_time:
                first_time = False

            if grps[iter] == [] and len(grps) != iter+1 and (not (any([len(x)>0 for x in out_exts]) and iter+1 > iter_max)):
                iter += 1
                if len(grps[iter]) > iter_item_cap:
                    logging.error("Recursion item cap hit not splitting {0} reactions at iter {1} with {2} items".format(len(template_rxn_map[parent.label]),iter,len(grps[iter])))
                    iter -= 1
                    gave_up_split = True

            elif grps[iter] == [] and len(grps) != iter+1 and (any([len(x)>0 for x in out_exts]) and iter+1 > iter_max):
                logging.error("iter_max achieved terminating early")

        out = []
        # compile all of the valid extensions together
        # may be some duplicates here, but I don't think it's currently worth identifying them
        for x in out_exts:
            out.extend(x)

        return out, gave_up_split

    def extend_node(self, parent, template_rxn_map, obj=None, T=1000.0, iter_max=np.inf, iter_item_cap=np.inf):
        """
        Constructs an extension to the group parent based on evaluation
        of the objective function obj
        """
        exts, gave_up_split = self.get_extension_edge(parent, template_rxn_map, obj=obj, T=T, iter_max=iter_max, iter_item_cap=iter_item_cap)

        if exts == [] and not gave_up_split:  # should only occur when all reactions at this node are identical
            rs = template_rxn_map[parent.label]
            for q, rxn in enumerate(rs):
                for j in range(q):
                    if not same_species_lists(rxn.reactants, rs[j].reactants, generate_initial_map=True, save_order=self.save_order):
                        for p, atm in enumerate(parent.item.atoms):
                            if atm.reg_dim_atm[0] != atm.reg_dim_atm[1]:
                                logging.error('atom violation')
                                logging.error(atm.reg_dim_atm)
                                logging.error(parent.label)
                                logging.error('Regularization dimension suggest this node can be expanded, '
                                              'but extension generation has failed')
                            if atm.reg_dim_u[0] != atm.reg_dim_u[1]:
                                logging.error('radical violation')
                                logging.error(atm.reg_dim_u)
                                logging.error(parent.label)
                                logging.error('Regularization dimension suggest this node can be expanded, '
                                              'but extension generation has failed')
                        for p, bd in enumerate(parent.item.get_all_edges()):
                            if bd.reg_dim[0] != bd.reg_dim[1]:
                                logging.error('bond violation')
                                logging.error(bd.order)
                                logging.error(bd.reg_dim)
                                logging.error(parent.label)
                                logging.error('Regularization dimension suggest this node can be expanded, '
                                              'but extension generation has failed')

                        logging.error('split violation')
                        logging.error('parent')
                        logging.error(parent.item.to_adjacency_list())
                        for c, atm in enumerate(parent.item.atoms):
                            logging.error(c)
                            logging.error(atm.reg_dim_atm)
                            logging.error(atm.reg_dim_u)
                        logging.error("bonds:")
                        for bd in parent.item.get_all_edges():
                            ind1 = parent.item.atoms.index(bd.vertex1)
                            ind2 = parent.item.atoms.index(bd.vertex2)
                            logging.error(((ind1, ind2), bd.order, bd.reg_dim))
                        for rxn in rs:
                            logging.error(str(rxn))
                            for react in rxn.reactants:
                                logging.error(react.label)
                                logging.error(react.to_adjacency_list())
                            for prod in rxn.products:
                                logging.error(prod.label)
                                logging.error(prod.to_adjacency_list())
                        logging.error("Clearing Regularization Dimensions and Reattempting")  # this usually happens when node expansion breaks some symmetry
                        parent.item.clear_reg_dims()  # this almost always solves the problem
                        return True
            return False

        if gave_up_split:
            return False

        vals = []
        for grp, grpc, name, typ, einds in exts:
            val, boo = self.eval_ext(parent, grp, name, template_rxn_map, obj, T)
            vals.append(val)

        min_val = min(vals)

        min_ind = vals.index(min_val)

        ext = exts[min_ind]

        extname = ext[2]

        if ext[3] == 'atomExt':
            ext[0].atoms[ext[4][0]].reg_dim_atm = [ext[0].atoms[ext[4][0]].atomtype, ext[0].atoms[ext[4][0]].atomtype]
        elif ext[3] == 'elExt':
            ext[0].atoms[ext[4][0]].reg_dim_u = [ext[0].atoms[ext[4][0]].radical_electrons,
                                                 ext[0].atoms[ext[4][0]].radical_electrons]

        self.add_entry(parent, ext[0], extname)

        complement = not ext[1] is None

        if complement:
            frags = extname.split('_')
            frags[-1] = 'N-' + frags[-1]
            cextname = ''
            for k in frags:
                cextname += k
                cextname += '_'
            cextname = cextname[:-1]

            self.add_entry(parent, ext[1], cextname)

        rxns = template_rxn_map[parent.label]

        comp_entries = []
        new_entries = []

        for i, entry in enumerate(template_rxn_map[parent.label]):
            if self.reaction_matches(entry,ext[0]):
                new_entries.append(entry)
            elif ext[1] is None or self.reaction_matches(entry,ext[1]):
                comp_entries.append(entry)
            else:
                logging.error("Reaction matched neither the new group or its complement")
                for p, atm in enumerate(parent.item.atoms):
                    if atm.reg_dim_atm[0] != atm.reg_dim_atm[1]:
                        logging.error('atom violation')
                        logging.error(atm.reg_dim_atm)
                        logging.error(parent.label)
                        logging.error('Regularization dimension suggest this node can be expanded, '
                                      'but extension generation has failed')
                    if atm.reg_dim_u[0] != atm.reg_dim_u[1]:
                        logging.error('radical violation')
                        logging.error(atm.reg_dim_u)
                        logging.error(parent.label)
                        logging.error('Regularization dimension suggest this node can be expanded, '
                                      'but extension generation has failed')
                for p, bd in enumerate(parent.item.get_all_edges()):
                    if bd.reg_dim[0] != bd.reg_dim[1]:
                        logging.error('bond violation')
                        logging.error(bd.order)
                        logging.error(bd.reg_dim)
                        logging.error(parent.label)
                        logging.error('Regularization dimension suggest this node can be expanded, '
                                      'but extension generation has failed')

                logging.error('split violation')
                logging.error("proposed")
                logging.error(ext[0].to_adjacency_list())
                if ext[1]:
                    logging.error("proposed complement")
                    logging.error(ext[1].to_adjacency_list())
                logging.error('parent')
                logging.error(parent.item.to_adjacency_list())
                for c, atm in enumerate(parent.item.atoms):
                    logging.error(c)
                    logging.error(atm.reg_dim_atm)
                    logging.error(atm.reg_dim_u)
                logging.error("bonds:")
                for bd in parent.item.get_all_edges():
                    ind1 = parent.item.atoms.index(bd.vertex1)
                    ind2 = parent.item.atoms.index(bd.vertex2)
                    logging.error(((ind1, ind2), bd.order, bd.reg_dim))
                rs = template_rxn_map[parent.label]
                for rxn in rs:
                    logging.error(str(rxn))
                    for react in rxn.reactants:
                        logging.error(react.label)
                        logging.error(react.to_adjacency_list())
                    for prod in rxn.products:
                        logging.error(prod.label)
                        logging.error(prod.to_adjacency_list())
                raise ValueError

        template_rxn_map[extname] = new_entries

        if complement:
            template_rxn_map[parent.label] = []
            template_rxn_map[cextname] = comp_entries
        else:
            template_rxn_map[parent.label] = comp_entries
        return True

    def generate_tree(self, rxns=None, obj=None, thermo_database=None, T=1000.0, nprocs=1, min_splitable_entry_num=2,
                      min_rxns_to_spawn=20, max_batch_size=800, outlier_fraction=0.02, stratum_num=8,
                      new_fraction_threshold_to_reopt_node=0.25, extension_iter_max=np.inf, extension_iter_item_cap=np.inf):
        """
        Generate a tree by greedy optimization based on the objective function obj
        the optimization is done by iterating through every group and if the group has
        more than one training reaction associated with it a set of potential more specific extensions
        are generated and the extension that optimizing the objective function combination is chosen
        and the iteration starts over at the beginning

        additionally the tree structure is simplified on the fly by removing groups that have no kinetics data
        associated if their parent has no kinetics data associated and they either have only one child or
        have two children one of which has no kinetics data and no children
        (its parent becomes the parent of its only relevant child node)

        Args:
            rxns: List of reactions to generate tree from (if None pull the whole training set)
            obj: Object to expand tree from (if None uses top node)
            thermo_database: Thermodynamic database used for reversing training reactions
            T: Temperature the tree is optimized for
            nprocs: Number of process for parallel tree generation
            min_splitable_entry_num: the minimum number of splitable reactions at a node in order to spawn
                a new process solving that node
            min_rxns_to_spawn: the minimum number of reactions at a node to spawn a new process solving that node
            max_batch_size: the maximum number of reactions allowed in a batch, most batches will be this size
                the last will be smaller, if the # of reactions < max_batch_size the cascade algorithm is not used
            outlier_fraction: Fraction of reactions that are fastest/slowest and will be automatically included
                in the first batch
            stratum_num: Number of strata used in stratified sampling scheme
            max_rxns_to_reopt_node: Nodes with more matching reactions than this will not be pruned
        """
        if rxns is None:
            rxns = self.get_training_set(thermo_database=thermo_database, remove_degeneracy=True, estimate_thermo=True,
                                         fix_labels=True, get_reverse=True, rxns_with_kinetics_only=True)

        if len(rxns) <= max_batch_size:
            template_rxn_map = self.get_reaction_matches(rxns=rxns, thermo_database=thermo_database, remove_degeneracy=True,
                                                         fix_labels=True, exact_matches_only=True, get_reverse=True)
            self.make_tree_nodes(template_rxn_map=template_rxn_map, obj=obj, T=T, nprocs=nprocs - 1, depth=0,
                                 min_splitable_entry_num=min_splitable_entry_num, min_rxns_to_spawn=min_rxns_to_spawn,extension_iter_max=extension_iter_max)
        else:
            def rxnkey(rxn):
                c = 0
                for react in rxn.reactants:
                    c += len(react.molecule[0].atoms)
                return c
            rxnsorted = sorted(rxns,key=rxnkey)
            batches = [rxnsorted[i * max_batch_size:(i + 1) * max_batch_size] for i in range((len(rxnsorted) + max_batch_size - 1) // max_batch_size )]
            for i, batch in enumerate(batches):
                if i == 0:
                    rxns = batch
                else:
                    rxns += batch
                    logging.error("pruning tree")
                    self.prune_tree(rxns, batch, thermo_database=thermo_database, new_fraction_threshold_to_reopt_node=new_fraction_threshold_to_reopt_node)
                    logging.error("pruned tree down to {} nodes".format(len(list(self.groups.entries))))
                logging.error("getting reaction matches")
                template_rxn_map = self.get_reaction_matches(rxns=rxns, thermo_database=thermo_database, fix_labels=True,
                                                             exact_matches_only=True, get_reverse=True)
                logging.error("building tree with {} rxns".format(len(rxns)))
                self.make_tree_nodes(template_rxn_map=template_rxn_map, obj=obj, T=T, nprocs=nprocs - 1, depth=0,
                                     min_splitable_entry_num=min_splitable_entry_num, min_rxns_to_spawn=min_rxns_to_spawn, extension_iter_max=extension_iter_max,
                                     extension_iter_item_cap=extension_iter_item_cap)
                logging.error("built tree with {} nodes".format(len(list(self.groups.entries))))

            self.auto_generated = True

    def get_rxn_batches(self, rxns, T=1000.0, max_batch_size=800, outlier_fraction=0.02, stratum_num=8):
        """
        Breaks reactions into batches based on a modified stratified sampling scheme
        Effectively:
        The top and bottom outlier_fraction of all reactions are always included in the first batch
        The remaining reactions are ordered by the rate coefficients at T
        The list of reactions is then split into stratum_num similarly sized intervals
        batches sample equally from each interval, but randomly within each interval
        until they reach max_batch_size reactions
        A list of lists of reactions containing the batches is returned
        """
        ks = np.array([rxn.kinetics.get_rate_coefficient(T=T) for rxn in rxns])
        inds = np.argsort(ks)
        outlier_num = int(outlier_fraction * len(ks) / 2)
        if outlier_num == 0:
            lowouts = []
            highouts = []
        else:
            lowouts = inds[:outlier_num].tolist()
            highouts = inds[-outlier_num:].tolist()
            inds = inds[outlier_num:-outlier_num]
        interval_length = int(len(inds) / stratum_num)
        strata = []
        for i in range(stratum_num):
            if i == 0:
                temp = inds[:interval_length].tolist()
                random.shuffle(temp)
                strata.append(temp)
            elif i == stratum_num - 1:
                temp = inds[interval_length * i:].tolist()
                random.shuffle(temp)
                strata.append(temp)
            else:
                temp = inds[interval_length * i:interval_length * (i + 1)].tolist()
                random.shuffle(temp)
                strata.append(temp)

        first_batch_strata_num = max_batch_size - outlier_num
        batches = [highouts + lowouts]
        bind = 0
        while any([len(stratum) != 0 for stratum in strata]):
            for stratum in strata:
                if stratum != []:
                    batches[bind].append(stratum.pop())
                    if len(batches[bind]) >= max_batch_size:
                        bind += 1
                        batches.append([])

        rxns = np.array(rxns)
        batches = [rxns[inds].tolist() for inds in batches if len(inds) > 0]

        return batches

    def prune_tree(self, rxns, newrxns, thermo_database=None, new_fraction_threshold_to_reopt_node=0.25, fix_labels=True,
                   exact_matches_only=True, get_reverse=True):
        """
        Remove nodes that have less than maxRxnToReoptNode reactions that match
        and clear the regularization dimensions of their parent
        This is used to remove smaller easier to optimize and more likely to change nodes
        before adding a new batch in cascade model generation
        """
        template_rxn_map = self.get_reaction_matches(rxns=rxns, thermo_database=thermo_database, fix_labels=fix_labels,
                                                     exact_matches_only=False, get_reverse=get_reverse)
        for key, item in template_rxn_map.items():
            entry = self.groups.entries[key]
            parent = entry.parent
            Nrxns = float(len(template_rxn_map[entry.label]))
            if parent and Nrxns > 0:
                Nrxnsnew = float(len(set(template_rxn_map[entry.label]) & set(newrxns)))
                if Nrxnsnew/Nrxns > new_fraction_threshold_to_reopt_node:
                    parent.children.remove(entry)
                    del self.groups.entries[key]
                else:
                    entry.item.clear_reg_dims()

    def make_tree_nodes(self, template_rxn_map=None, obj=None, T=1000.0, nprocs=0, depth=0, min_splitable_entry_num=2,
                        min_rxns_to_spawn=20, extension_iter_max=np.inf, extension_iter_item_cap=np.inf):

        if depth > 0:
            root = self.groups.entries[list(template_rxn_map.keys())[0]]
        else:
            for entry in self.groups.entries.values():  # find the root entry for this branch
                if entry.index != -1:
                    root = entry
                    break
            while root.parent is not None:
                root = root.parent

        if depth == 0: #may start with reactions at many nodes first iteration
            psize = 0.0
            entries = [root]
            while entries != []:
                entry = entries[-1]
                if entry.children:
                    psize += float(len(template_rxn_map[entry.label]))
                    entries.extend(entry.children)
                    entries.remove(entry)
                else:
                    psize += float(len(template_rxn_map[entry.label]))
                    entries.remove(entry)
        else:
            psize = float(len(template_rxn_map[root.label]))

        logging.error(psize)
        mult_completed_nodes = []  # nodes containing multiple identical training reactions
        boo = True  # if the for loop doesn't break becomes false and the while loop terminates
        active_procs = []
        active_conns = []
        active_proc_num = []
        proc_names = []
        free_procs = nprocs
        extra_entries = []

        while boo:
            remove_inds = []
            for k, p in enumerate(active_procs):  # check if any processes have finished
                if active_conns[k].poll():
                    new_entries = self._absorb_process(p, active_conns[k], proc_names[k])
                    extra_entries += new_entries
                    remove_inds.append(k)

            remove_inds.reverse()
            for ind in remove_inds:  # remove finished process objects
                free_procs += active_proc_num[ind]
                del active_proc_num[ind]
                del active_procs[ind]
                del active_conns[ind]
                del proc_names[ind]

            splitable_entry_num = 0
            for label, items in template_rxn_map.items():  # figure out how many splitable objects there are
                entry = self.groups.entries[label]
                if len(items) > 1 and entry not in mult_completed_nodes:
                    splitable_entry_num += 1

            for label in list(template_rxn_map.keys()):
                entry = self.groups.entries[label]
                if not isinstance(entry.item, Group):  # skip logic nodes
                    continue
                if len(template_rxn_map[label]) == 0:
                    continue
                if entry.index != -1 and len(template_rxn_map[entry.label]) > 1 and entry not in mult_completed_nodes:
                    if entry.parent and (free_procs > 0 and splitable_entry_num > min_splitable_entry_num and
                            len(template_rxn_map[entry.label]) > min_rxns_to_spawn):
                        procs_out = int(len(template_rxn_map[entry.label]) / psize * free_procs)
                        free_procs -= procs_out
                        assert free_procs >= 0
                        conn, p, name = _spawn_tree_process(family=self,
                                                            template_rxn_map={entry.label: template_rxn_map[entry.label]},
                                                            obj=obj, T=T, nprocs=procs_out - 1, depth=depth,
                                                            min_splitable_entry_num=min_splitable_entry_num,
                                                            min_rxns_to_spawn=min_rxns_to_spawn,extension_iter_max=extension_iter_max,
                                                            extension_iter_item_cap=extension_iter_item_cap)
                        active_procs.append(p)
                        active_conns.append(conn)
                        proc_names.append(name)
                        active_proc_num.append(procs_out)
                        L = entry.label
                        self.groups.entries[L].parent.children.remove(entry)
                        del template_rxn_map[L]  # prevents this process from recreating work done by another process
                        del self.groups.entries[L]

                        splitable_entry_num -= 1
                        continue
                    boo2 = self.extend_node(entry, template_rxn_map, obj, T, iter_max=extension_iter_max, iter_item_cap=extension_iter_item_cap)
                    if boo2:  # extended node so restart while loop
                        break
                    else:  # no extensions could be generated since all reactions were identical
                        mult_completed_nodes.append(entry)
            else:
                if len(active_procs) == 0:
                    boo = False

            # fix indicies
            iters = 0
            for entry in self.groups.entries.values():
                if entry.index != -1:
                    entry.index = iters
                    iters += 1

        # add the entries generated on other processors
        index = max([ent.index for ent in self.groups.entries.values()]) + 1
        for item in extra_entries:
            if item.label in self.groups.entries.keys():
                continue
            item.index = index
            index += 1
            self.groups.entries[item.label] = item

        for label, entry in self.groups.entries.items():
            if entry.index != -1 and entry.parent is None and entry.label != root.label:
                pname = "_".join(label.split('_')[:-1])
                while pname not in self.groups.entries.keys():
                    pname = "_".join(label.split('_')[:-1])
                entry.parent = self.groups.entries[pname]
                entry.parent.children.append(entry)

        return

    def _absorb_process(self, p, conn, name):
        try:
            grps = conn.recv()
            p.terminate()
        except Exception as e:
            logging.error('failed to absorb process {}'.format(name))
            raise e
        return grps

    def make_bm_rules_from_template_rxn_map(self, template_rxn_map, nprocs=1, Tref=1000.0, fmax=1.0e5):

        rule_keys = self.rules.entries.keys()
        for entry in self.groups.entries.values():
            if entry.label not in rule_keys:
                self.rules.entries[entry.label] = []

        index = max([e.index for e in self.rules.get_entries()] or [0]) + 1

        entries = list(self.groups.entries.values())
        rxnlists = [(template_rxn_map[entry.label], entry.label)
                    if entry.label in template_rxn_map.keys() else [] for entry in entries]
        inputs = [(self.forward_recipe.actions, rxns, Tref, fmax, label, [r.rank for r in rxns])
                           for rxns, label in rxnlists]

        inds = np.arange(len(inputs))
        np.random.shuffle(inds)  # want to parallelize in random order
        inds = inds.tolist()
        revinds = [inds.index(x) for x in np.arange(len(inputs))]

        if nprocs > 1:
            pool = mp.Pool(nprocs)
            kinetics_list = np.array(pool.map(_make_rule, list(inputs[i] for i in inds)))
        else:
            kinetics_list = np.array(list(map(_make_rule, list(inputs[i] for i in inds))))

        kinetics_list = kinetics_list[revinds]  # fix order

        for i, kinetics in enumerate(kinetics_list):
            if isinstance(kinetics, Marcus):
                entry = entries[i]
                st = "Marcus rule fitted to {0} training reactions at node {1}".format(len(rxnlists[i][0]), entry.label)
                new_entry = Entry(
                    index=index,
                    label=entry.label,
                    item=self.forward_template,
                    data=kinetics,
                    rank=11,
                    reference=None,
                    short_desc=st,
                    long_desc=st,
                )
                new_entry.data.comment = st

                self.rules.entries[entry.label].append(new_entry)

                index += 1
            elif kinetics is not None:
                entry = entries[i]
                std = kinetics.uncertainty.get_expected_log_uncertainty() / 0.398  # expected uncertainty is std * 0.398
                st = "BM rule fitted to {0} training reactions at node {1}".format(len(rxnlists[i][0]), entry.label)
                st += "\nTotal Standard Deviation in ln(k): {0}".format(std)
                new_entry = Entry(
                    index=index,
                    label=entry.label,
                    item=self.forward_template,
                    data=kinetics,
                    rank=11,
                    reference=None,
                    short_desc=st,
                    long_desc=st,
                )
                new_entry.data.comment = st

                self.rules.entries[entry.label].append(new_entry)

                index += 1

        for label,entry in self.rules.entries.items(): #pull solute data from further up the tree as needed
            if len(entry) == 0:
                continue
            entry = entry[0]
            if not entry.data.solute:
                ent = self.groups.entries[label]
                while not self.rules.entries[ent.label][0].data.solute and ent.parent:
                    ent = ent.parent
                entry.data.solute = self.rules.entries[ent.label][0].data.solute

    def cross_validate(self, folds=5, template_rxn_map=None, test_rxn_inds=None, T=1000.0, iters=0, random_state=1):
        """
        Perform K-fold cross validation on an automatically generated tree at temperature T
        after finding an appropriate node for kinetics estimation it will move up the tree
        iters times.
        Returns a dictionary mapping {rxn:Ln(k_Est/k_Train)}
        """

        if template_rxn_map is None:
            template_rxn_map = self.get_reaction_matches(remove_degeneracy=True, get_reverse=True, fix_labels=True)

        rxns = np.array(template_rxn_map['Root'])

        if test_rxn_inds is None:
            if folds == 0:
                folds = len(rxns)

            kf = KFold(folds, shuffle=True, random_state=random_state)
            kfsplits = kf.split(rxns)
        else:
            kfsplits = [([0, ], [0, ])]

        errors = {}
        uncertainties = {}

        for train_index, test_index in kfsplits:

            if test_rxn_inds is None:
                rxns_test = rxns[test_index]
            else:
                rxns_test = rxns[test_rxn_inds]

            for rxn in rxns_test:

                krxn = rxn.kinetics.get_rate_coefficient(T)

                entry = self.get_root_template()[0]

                boo = True
                while boo:  # find the entry it matches
                    for child in entry.children:
                        rs = template_rxn_map[child.label]
                        if rxn in rs:
                            entry = child
                            break
                    else:
                        boo = False

                while entry.parent and len(set(template_rxn_map[entry.label]) - set(rxns_test)) <= 1:
                    if entry.parent:
                        entry = entry.parent

                for q in range(iters):
                    if entry.parent:
                        entry = entry.parent

                boo = True

                while boo:
                    if entry.parent is None:
                        break
                    kin = self.rules.entries[entry.label][0].data
                    kinparent = self.rules.entries[entry.parent.label][0].data
                    err_parent = abs(kinparent.uncertainty.data_mean + kinparent.uncertainty.mu - kin.uncertainty.data_mean) + np.sqrt(2.0*kinparent.uncertainty.var/np.pi)
                    err_entry = abs(kin.uncertainty.mu) + np.sqrt(2.0*kin.uncertainty.var/np.pi)
                    if err_entry <= err_parent:
                        break
                    else:
                        entry = entry.parent

                uncertainties[rxn] = self.rules.entries[entry.label][0].data.uncertainty


                L = list(set(template_rxn_map[entry.label]) - set(rxns_test))

                if L != []:
                    if isinstance(L[0].kinetics, Arrhenius):
                        kinetics = ArrheniusBM().fit_to_reactions(L, recipe=self.forward_recipe.actions)
                        if kinetics.E0.value_si < 0.0 or len(L) == 1:
                            kinetics = average_kinetics([r.kinetics for r in L])
                        else:
                            kinetics = kinetics.to_arrhenius(rxn.get_enthalpy_of_reaction(298.))
                    else:
                        kinetics = ArrheniusChargeTransferBM().fit_to_reactions(L, recipe=self.forward_recipe.actions)
                        if kinetics.E0.value_si < 0.0 or len(L) == 1:
                            kinetics = average_kinetics([r.kinetics for r in L])
                        else:
                            kinetics = kinetics.to_arrhenius_charge_transfer(rxn.get_enthalpy_of_reaction(298.))

                    k = kinetics.get_rate_coefficient(T)
                    errors[rxn] = np.log(k / krxn)
                else:
                    raise ValueError('only one piece of kinetics information in the tree?')

        return errors, uncertainties

    def cross_validate_old(self, folds=5, T=1000.0, random_state=1, estimator='rate rules', thermo_database=None, get_reverse=False, uncertainties=True):
        """
        Perform K-fold cross validation on an automatically generated tree at temperature T
        Returns a dictionary mapping {rxn:Ln(k_Est/k_Train)}
        """
        errors = {}
        uncs = {}

        kpu = KineticParameterUncertainty()
        rxns = np.array(self.get_training_set(remove_degeneracy=True,get_reverse=get_reverse))

        if folds == 0:
            folds = len(rxns)

        kf = KFold(folds, shuffle=True, random_state=random_state)

        if thermo_database is None:
            from rmgpy.data.rmg import get_db
            tdb = get_db('thermo')
        else:
            tdb = thermo_database

        for train_index, test_index in kf.split(rxns):

            self.rules.entries = {}  # clear rules each iteration
            if get_reverse:
                train_index = [x for x in train_index if x<len(rxns)/2]
            self.add_rules_from_training(train_indices=train_index, thermo_database=tdb)
            self.fill_rules_by_averaging_up()
            rxns_test = rxns[test_index]

            for rxn in rxns_test:

                krxn = rxn.kinetics.get_rate_coefficient(T)

                template_labels = self.get_reaction_template_labels(rxn)
                template = self.retrieve_template(template_labels)
                if estimator == 'rate rules':
                    kinetics, entry = self.estimate_kinetics_using_rate_rules(template, degeneracy=1)
                else:
                    raise ValueError('{0} is not a valid value for input `estimator`'.format(estimator))

                k = kinetics.get_rate_coefficient(T)

                errors[rxn] = np.log(k / krxn)
                if uncertainties:
                    testrxn = deepcopy(rxn)
                    testrxn.kinetics = kinetics
                    boo,source = self.extract_source_from_comments(testrxn)
                    sdict = {"Rate Rules":source}
                    uncs[rxn] = kpu.get_uncertainty_value(sdict)

        if uncertainties:
            return errors, uncs
        else:
            return errors

    def simple_regularization(self, node, template_rxn_map, test=True):
        """
        Simplest regularization algorithm
        All nodes are made as specific as their descendant reactions
        Training reactions are assumed to not generalize
        For example if an particular atom at a node is Oxygen for all of its
        descendent reactions a reaction where it is Sulfur will never hit that node
        unless it is the top node even if the tree did not split on the identity
        of that atom

        The test option to this function determines whether or not the reactions
        under a node match the extended group before adding an extension.
        If the test fails the extension is skipped.

        In general test=True is needed if the cascade algorithm was used
        to generate the tree and test=False is ok if the cascade algorithm
        wasn't used.
        """

        for child in node.children:
            self.simple_regularization(child, template_rxn_map)

        grp = node.item
        rxns = template_rxn_map[node.label]

        R = ['H', 'C', 'N', 'O', 'Si', 'S', 'Cl', 'F', 'Br']  # set of possible R elements/atoms
        R = [ATOMTYPES[x] for x in R]

        RnH = R[:]
        RnH.remove(ATOMTYPES['H'])

        Run = [0, 1, 2, 3]

        atm_dict = {'R': R, 'R!H': RnH}

        if isinstance(node.item, Group):
            indistinguishable = []
            for i, atm1 in enumerate(grp.atoms):

                skip = False
                if node.children == []:  # if the atoms or bonds are graphically indistinguishable don't regularize
                    bdpairs = {(atm, tuple(bd.order)) for atm, bd in atm1.bonds.items()}
                    for atm2 in grp.atoms:
                        if atm1 is not atm2 and atm1.atomtype == atm2.atomtype and len(atm1.bonds) == len(atm2.bonds):
                            bdpairs2 = {(atm, tuple(bd.order)) for atm, bd in atm2.bonds.items()}
                            if bdpairs == bdpairs2:
                                skip = True
                                indistinguishable.append(i)

                if not skip and atm1.reg_dim_atm[1] != [] and set(atm1.reg_dim_atm[1]) != set(atm1.atomtype):
                    atyp = atm1.atomtype
                    if len(atyp) == 1 and atyp[0] in R:
                        pass
                    else:
                        if len(atyp) == 1 and atyp[0].label in atm_dict.keys():
                            atyp = atm_dict[atyp[0].label]

                        vals = list(set(atyp) & set(atm1.reg_dim_atm[1]))
                        assert vals != [], 'cannot regularize to empty'
                        if all([set(child.item.atoms[i].atomtype) <= set(vals) for child in node.children]):
                            if not test:
                                atm1.atomtype = vals
                            else:
                                oldvals = atm1.atomtype
                                atm1.atomtype = vals
                                if not self.rxns_match_node(node, rxns):
                                    atm1.atomtype = oldvals

                if not skip and atm1.reg_dim_u[1] != [] and set(atm1.reg_dim_u[1]) != set(atm1.radical_electrons):
                    if len(atm1.radical_electrons) == 1:
                        pass
                    else:
                        relist = atm1.radical_electrons
                        if relist == []:
                            relist = Run
                        vals = list(set(relist) & set(atm1.reg_dim_u[1]))
                        assert vals != [], 'cannot regularize to empty'

                        if all([set(child.item.atoms[i].radical_electrons) <= set(vals)
                                if child.item.atoms[i].radical_electrons != [] else False for child in node.children]):
                            if not test:
                                atm1.radical_electrons = vals
                            else:
                                oldvals = atm1.radical_electrons
                                atm1.radical_electrons = vals
                                if not self.rxns_match_node(node, rxns):
                                    atm1.radical_electrons = oldvals

                if (not skip and atm1.reg_dim_r[1] != [] and
                        ('inRing' not in atm1.props.keys() or atm1.reg_dim_r[1][0] != atm1.props['inRing'])):
                    if 'inRing' not in atm1.props.keys():
                        if (all(['inRing' in child.item.atoms[i].props.keys() for child in node.children]) and
                                all([child.item.atoms[i].props['inRing'] == atm1.reg_dim_r[1] for child in node.children])):
                            if not test:
                                atm1.props['inRing'] = atm1.reg_dim_r[1][0]
                            else:
                                if 'inRing' in atm1.props.keys():
                                    oldvals = atm1.props['inRing']
                                else:
                                    oldvals = None
                                atm1.props['inRing'] = atm1.reg_dim_r[1][0]
                                if not self.rxns_match_node(node, rxns):
                                    if oldvals:
                                        atm1.props['inRing'] = oldvals
                                    else:
                                        del atm1.props['inRing']
                if not skip:
                    for j, atm2 in enumerate(grp.atoms[:i]):
                        if j in indistinguishable:  # skip graphically indistinguishable atoms
                            continue
                        if grp.has_bond(atm1, atm2):
                            bd = grp.get_bond(atm1, atm2)
                            if len(bd.order) == 1:
                                pass
                            else:
                                vals = list(set(bd.order) & set(bd.reg_dim[1]))
                                if vals != [] and all([set(child.item.get_bond(child.item.atoms[i], child.item.atoms[j]).order) <= set(vals) for child in node.children]):
                                    if not test:
                                        bd.order = vals
                                    else:
                                        oldvals = bd.order
                                        bd.order = vals
                                        if not self.rxns_match_node(node, rxns):
                                            bd.order = oldvals

    def regularize(self, regularization=simple_regularization, keep_root=True, thermo_database=None,
                   template_rxn_map=None, rxns=None):
        """
        Regularizes the tree according to the regularization function regularization
        """
        if template_rxn_map is None:
            if rxns is None:
                template_rxn_map = self.get_reaction_matches(thermo_database=thermo_database, remove_degeneracy=True,
                                                             get_reverse=True, exact_matches_only=False, fix_labels=True, rxns_with_kinetics_only=False)
            else:
                template_rxn_map = self.get_reaction_matches(rxns=rxns, thermo_database=thermo_database,
                                                             remove_degeneracy=True, get_reverse=True, exact_matches_only=False,
                                                             fix_labels=True)

        if keep_root:
            for child in self.get_root_template()[0].children:  # don't regularize the root
                regularization(self, child, template_rxn_map)
        else:
            regularization(self, self.get_root_template()[0], template_rxn_map)

    def check_tree(self, entry=None):
        if entry is None:
            entry = self.get_root_template()[0]
        for child in entry.children:
            if not child.item.is_subgraph_isomorphic(entry.item, generate_initial_map=True, save_order=True):
                logging.error('child: ')
                logging.error(child.label)
                logging.error(child.item.to_adjacency_list())
                logging.error('parent: ')
                logging.error(entry.label)
                logging.error(entry.item.to_adjacency_list())
                raise ValueError('Child not subgraph isomorphic to parent')
            self.check_tree(child)
        for entry in self.groups.entries.values():
            if entry.index == -1:
                continue
            parent = entry
            while parent.parent is not None:
                parent = parent.parent
            assert parent.label == 'Root', parent.label

    def make_tree(self, obj=None, regularization=simple_regularization, thermo_database=None, T=1000.0):
        """
        generates tree structure and then generates rules for the tree
        """
        self.generate_tree(obj=obj, thermo_database=thermo_database, T=T)
        self.regularize(regularization=regularization)
        template_rxn_map = self.get_reaction_matches(thermo_database=thermo_database,
                                                     remove_degeneracy=True, get_reverse=True)
        self.make_bm_rules_from_template_rxn_map(template_rxn_map)
        self.check_tree()

    def clean_tree_rules(self):
        self.rules.entries = OrderedDict()
        self.rules.entries['Root'] = []

    def clean_tree_groups(self):
        """
        clears groups and rules in the tree, generates an appropriate
        root group to start from and then reads training reactions
        Note this only works if a single top node (not a logic node)
        can be generated
        """
        # find the starting node
        grp = None

        rtmps = self.get_root_template()

        if not isinstance(rtmps[0].item, Group):
            raise ValueError('each tree top node must be a group not a logic node to prepare the tree automatically')

        for ent in rtmps:
            if grp is None:
                grp = ent.item
            else:
                if any([isinstance(x, list) for x in ent.item.get_all_labeled_atoms().values()]):
                    grp = grp.merge_groups(ent.item, keep_identical_labels=True)
                else:
                    grp = grp.merge_groups(ent.item)

        # clear everything
        self.groups.entries = {x.label: x for x in self.groups.entries.values() if x.index == -1}

        # add the starting node
        self.add_entry(None, grp, 'Root')
        self.groups.entries['Root'].index = 1
        self.groups.top = [self.groups.entries['Root']]
        self.forward_template.reactants = [self.groups.entries['Root']]

        return

    def clean_tree(self):
        self.clean_tree_rules()
        self.clean_tree_groups()

    def save_generated_tree(self, path=None):
        """
        clears the rules and saves the family to its
        current location in database
        """
        if path is None:
            path = settings['database.directory']
            path = os.path.join(path, 'kinetics', 'families', self.label)

        self.save(path)

    def get_training_set(self, thermo_database=None, remove_degeneracy=False, estimate_thermo=True, fix_labels=False,
                         get_reverse=False, rxns_with_kinetics_only=False):
        """
        retrieves all reactions in the training set, assigns thermo to the species objects
        reverses reactions as necessary so that all reactions are in the forward direction
        and returns the resulting list of reactions in the forward direction with thermo
        assigned
        """

        def get_label_fixed_mol(mol, root_labels):
            nmol = mol.copy(deep=True)
            for atm in nmol.atoms:
                if atm.label not in root_labels:
                    atm.label = ''
            return nmol

        def fix_labels_mol(mol, root_labels):
            for atm in mol.atoms:
                if atm.label not in root_labels:
                    atm.label = ''

        def get_reactant_thermo(reactant,metal):
            """
            Save the label of reactant and reapply them after thermo estimation to avoid deepcopying
            """
            label_dict = reactant.molecule[0].get_all_labeled_atoms()
            mol = reactant.molecule[0]
            reactant.molecule[0].clear_labeled_atoms()
            reactant.generate_resonance_structures()
            if metal:
                thermo = tdb.get_thermo_data(reactant, metal_to_scale_to=metal)
            else:
                thermo = tdb.get_thermo_data(reactant)
            reactant.molecule = [mol]

            for key,atm in label_dict.items():
                if isinstance(atm,list):
                    for a in atm:
                        a.label = key
                else:
                    atm.label = key
            return thermo

        if self.own_reverse and get_reverse:
            rev_rxns = []
            rkeys = list(self.reverse_map.keys())
            reverse_map = self.reverse_map

        if estimate_thermo:
            if thermo_database is None:
                from rmgpy.data.rmg import get_db
                tdb = get_db('thermo')
            else:
                tdb = thermo_database

        try:
            dep = self.get_training_depository()
        except:
            logging.info('Could not find training depository in family {0}.'.format(self.label))
            logging.info('Must be because you turned off the training depository.')
            return

        rxns = deepcopy([i.item for i in dep.entries.values() if (not rxns_with_kinetics_only) or type(i.data) != KineticsModel])
        entries = deepcopy([i for i in dep.entries.values() if (not rxns_with_kinetics_only) or type(i.data) != KineticsModel])

        roots = [x.item for x in self.get_root_template()]
        root = None
        for r in roots:
            if root:
                root = root.merge_groups(r)
            else:
                root = deepcopy(r)

        root_labels = [x.label for x in root.atoms if x.label != '']
        root_label_set = set(root_labels)

        for i, entry in enumerate(entries):
            if estimate_thermo:
                # parse out the metal to scale to
                if entry.facet is None:
                    metal = entry.metal # could be None
                else:
                    metal = entry.metal + entry.facet
                for j, react in enumerate(entry.item.reactants):
                    if rxns[i].reactants[j].thermo is None:
                        rxns[i].reactants[j].thermo = get_reactant_thermo(react,metal)

                for j, react in enumerate(entry.item.products):
                    if rxns[i].products[j].thermo is None:
                        rxns[i].products[j].thermo = get_reactant_thermo(react,metal)
            rxns[i].kinetics = entry.data
            rxns[i].rank = entry.rank

            if remove_degeneracy and type(rxns[i].kinetics) != KineticsModel:  # adjust for degeneracy
                rxns[i].kinetics.A.value_si /= rxns[i].degeneracy

            mol = None
            for react in rxns[i].reactants:
                if fix_labels:
                    fix_labels_mol(react.molecule[0], root_labels)
                if mol:
                    mol = mol.merge(react.molecule[0])
                else:
                    mol = deepcopy(react.molecule[0])

            mol.update_atomtypes()

            if fix_labels:
                for prod in rxns[i].products:
                    fix_labels_mol(prod.molecule[0], root_labels)
                for atm in mol.atoms:
                    if atm.label not in root_labels:
                        atm.label = ''

            if (mol.is_subgraph_isomorphic(root, generate_initial_map=True) or
                    (not fix_labels and
                     get_label_fixed_mol(mol, root_labels).is_subgraph_isomorphic(root, generate_initial_map=True))):
                rxns[i].is_forward = True
                if self.own_reverse and get_reverse:
                    mol = None
                    for react in rxns[i].products:
                        if mol:
                            mol = mol.merge(react.molecule[0])
                        else:
                            mol = deepcopy(react.molecule[0])

                    if fix_labels:
                        mol_label_set = set([x.label for x in get_label_fixed_mol(mol, root_labels).atoms if x.label != ''])
                    else:
                        mol_label_set = set([x.label for x in mol.atoms if x.label != ''])

                    if mol_label_set == root_label_set and ((mol.is_subgraph_isomorphic(root, generate_initial_map=True) or
                            (not fix_labels and
                             get_label_fixed_mol(mol, root_labels).is_subgraph_isomorphic(root, generate_initial_map=True)))):
                        # try product structures
                        products = [Species(molecule=[get_label_fixed_mol(x.molecule[0], root_labels)], thermo=x.thermo)
                                    for x in rxns[i].products]
                    else:
                        products = self.apply_recipe([s.molecule[0] for s in rxns[i].reactants], forward=True)
                        products = [Species(molecule=[p]) for p in products]

                    prodmol = None
                    for react in rxns[i].products:
                        if prodmol:
                            prodmol = prodmol.merge(react.molecule[0])
                        else:
                            prodmol = deepcopy(react.molecule[0])

                    if not prodmol.is_subgraph_isomorphic(root, generate_initial_map=True):
                        mol = None
                        for react in products:
                            if mol:
                                mol = mol.merge(react.molecule[0])
                            else:
                                mol = deepcopy(react.molecule[0])
                        if not mol.is_subgraph_isomorphic(root, generate_initial_map=True):
                            for p in products:
                                for atm in p.molecule[0].atoms:
                                    if atm.label in rkeys:
                                        atm.label = reverse_map[atm.label]

                    reacts = [Species(molecule=[get_label_fixed_mol(x.molecule[0], root_labels)], thermo=x.thermo)
                              for x in rxns[i].reactants]
                    if type(rxns[i].kinetics) != KineticsModel:
                        if rxns[i].kinetics.solute:
                            rxns[i].kinetics.solute = to_soluteTSdata(rxns[i].kinetics.solute,reactants=rxns[i].reactants)
                        rrev = Reaction(reactants=products, products=reacts,
                                    kinetics=rxns[i].generate_reverse_rate_coefficient(), rank=rxns[i].rank)
                    rrev.is_forward = False

                    if estimate_thermo:
                        for rev_react in rrev.reactants:
                            if rev_react.thermo is None:
                                rev_react.thermo = get_reactant_thermo(rev_react,metal)

                    rev_rxns.append(rrev)

                continue
            else:
                if self.own_reverse:
                    logging.error("rxn")
                    logging.error(str(rxns[i]))
                    logging.error("root")
                    logging.error(root.to_adjacency_list())
                    logging.error("mol")
                    logging.error(mol.to_adjacency_list())
                    raise ValueError("couldn't match reaction")

                mol = None
                for react in rxns[i].products:
                    if mol:
                        mol = mol.merge(react.molecule[0])
                    else:
                        mol = deepcopy(react.molecule[0])

                mol.update_atomtypes()

                if (mol.is_subgraph_isomorphic(root, generate_initial_map=True) or
                        (not fix_labels and
                         get_label_fixed_mol(mol, root_labels).is_subgraph_isomorphic(root, generate_initial_map=True))):  # try product structures
                    products = [Species(molecule=[get_label_fixed_mol(x.molecule[0], root_labels)], thermo=x.thermo)
                                for x in rxns[i].products]
                else:
                    products = self.apply_recipe([s.molecule[0] for s in rxns[i].reactants], forward=True)
                    products = [Species(molecule=[p]) for p in products]

                rrev = Reaction(reactants=products, products=rxns[i].reactants,
                                kinetics=rxns[i].generate_reverse_rate_coefficient(), rank=rxns[i].rank)

                rrev.is_forward = False

                if estimate_thermo:
                    for rev_react in rrev.reactants:
                        if rev_react.thermo is None:
                            rev_react.thermo = get_reactant_thermo(rev_react,metal)
                rxns[i] = rrev

        if self.own_reverse and get_reverse:
            return rxns + rev_rxns
        else:
            return rxns

    def get_reaction_matches(self, rxns=None, thermo_database=None, remove_degeneracy=False, estimate_thermo=True,
                             fix_labels=False, exact_matches_only=False, get_reverse=False, rxns_with_kinetics_only=False):
        """
        returns a dictionary mapping for each entry in the tree:
        (entry.label,entry.item) : list of all training reactions (or the list given) that match that entry
        """
        if rxns is None:
            rxns = self.get_training_set(thermo_database=thermo_database, remove_degeneracy=remove_degeneracy,
                                         estimate_thermo=estimate_thermo, fix_labels=fix_labels,
                                         get_reverse=get_reverse,rxns_with_kinetics_only=rxns_with_kinetics_only)

        entries = self.groups.entries

        assert len(set(entries.keys())) == len(entries.keys()), 'there are duplicate indices in family.group.entries'

        rxn_lists = {entry.label: [] for entry in entries.values()}

        root = self.get_root_template()[0]

        for rxn in rxns:
            mol = None
            for r in rxn.reactants:
                if mol is None:
                    mol = deepcopy(r.molecule[0])
                else:
                    mol = mol.merge(r.molecule[0])
            try:
                flag = not self.is_entry_match(mol, root, resonance=True)
            except:
                flag = not self.is_entry_match(mol, root, resonance=False)

            if flag:
                logging.error(root.item.to_adjacency_list())
                logging.error(mol.to_adjacency_list())
                for r in rxn.reactants:
                    logging.error(r.molecule[0].to_adjacency_list())
                for r in rxn.products:
                    logging.error(r.molecule[0].to_adjacency_list())
                raise ValueError('reaction: {0} does not match root template in family {1}'.format(rxn, self.label))

            rxn_lists[root.label].append(rxn)

            entry = root

            while entry.children != []:
                for child in entry.children:
                    if self.is_entry_match(mol, child, resonance=False):
                        entry = child
                        rxn_lists[child.label].append(rxn)
                        break
                else:
                    break

        if exact_matches_only:
            new_lists = dict()
            for key, rs in rxn_lists.items():
                newrs = set(rs)
                for child in self.groups.entries[key].children:
                    newrs -= set(rxn_lists[child.label])
                new_lists[key] = list(newrs)
            rxn_lists = new_lists

        return rxn_lists

    def is_entry_match(self, mol, entry, resonance=True):
        """
        determines if the labeled molecule object of reactants matches the entry entry
        """
        if isinstance(entry.item, Group):
            if resonance:
                structs = mol.generate_resonance_structures()
            else:
                structs = [mol]
            return any([mol.is_subgraph_isomorphic(entry.item, generate_initial_map=True) for mol in structs])
        elif isinstance(entry.item, LogicOr):
            return any([self.is_entry_match(mol, self.groups.entries[c], resonance=resonance)
                        for c in entry.item.components])

    def rxns_match_node(self, node, rxns):
        for rxn in rxns:
            mol = None
            for r in rxn.reactants:
                if mol is None:
                    mol = deepcopy(r.molecule[0])
                else:
                    mol = mol.merge(r.molecule[0])

            if not self.is_entry_match(mol, node, resonance=False):
                return False

        return True

    def retrieve_original_entry(self, template_label):
        """
        Retrieves the original entry, be it a rule or training reaction, given
        the template label in the form 'group1;group2' or 'group1;group2;group3'

        Returns tuple in the form
        (RateRuleEntry, TrainingReactionEntry)

        Where the TrainingReactionEntry is only present if it comes from a training reaction
        """
        template_labels = template_label.split()[-1].split(';')
        template = self.retrieve_template(template_labels)
        rule = self.get_rate_rule(template)
        if 'From training reaction' in rule.data.comment:
            training_index = int(rule.data.comment.split()[3])
            training_depository = self.get_training_depository()
            return rule, training_depository.entries[training_index]
        else:
            return rule, None

    def get_sources_for_template(self, template):
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

        def assign_weights_to_entries(entry_nested_list, weighted_entries, n=1):
            """
            Assign weights to an average of average nested list. Where n is the
            number of values being averaged recursively.
            """
            n = len(entry_nested_list) * n
            for entry in entry_nested_list:
                if isinstance(entry, list):
                    assign_weights_to_entries(entry, weighted_entries, n)
                else:
                    weighted_entries.append((entry, 1 / float(n)))
            return weighted_entries

        kinetics, entry = self.estimate_kinetics_using_rate_rules(template)
        if entry:
            return [(entry, 1)], []  # Must be a rate rule
        else:
            # The template was estimated using an average or another node
            rules = []
            training = []

            lines = kinetics.comment.split('\n')

            # remove the Euclidean distance and family lines to help parser
            lines = [line for line in lines if not line.startswith('Euclid') and not line.startswith('family:')]

            # Discard the last line, unless it's the only line!
            # The last line is 'Estimated using ... for rate rule (originalTemplate)'
            # if from training reaction is in the first line append it to the end of the second line and skip the first line
            if 'Average of' not in kinetics.comment:
                if 'From training reaction' in lines[0]:
                    comment = lines[1]
                else:
                    comment = lines[0]
                if comment.startswith('Estimated using template'):
                    token_template_label = comment.split()[3][1:-1]
                    rule_entry, training_entry = self.retrieve_original_entry(token_template_label)
                    if training_entry:
                        training.append((rule_entry, training_entry, 1))  # Weight is 1
                    else:
                        rules.append((rule_entry, 1))
                else:
                    raise ValueError('Could not parse unexpected line found in kinetics comment: {}'.format(comment))
            else:
                comment = ' '.join(lines[:-1])
                # Clean up line for exec
                eval_comment_string = \
                    re.sub(r" \+ ", ",",                        # any remaining + signs
                    re.sub(r"Average of ", "",                  # average of averages
                    re.sub(r"Average of \[(?!Average)", "['",   # average of groups
                    re.sub(r"(\w|\))]", r"\1']",                # initial closing bracket
                    re.sub(r"(?<=[\w)]) \+ (?=Average)", "',",  # + sign between non-average and average
                    re.sub(r"(?<=]) \+ (?!Average)", ",'",      # + sign between average and non-average
                    re.sub(r"(?<!]) \+ (?!Average)", "','",     # + sign between non-averages
                    comment)))))))

                entry_nested_list = eval(eval_comment_string)

                weighted_entries = assign_weights_to_entries(entry_nested_list, [])

                rules = {}
                training = {}

                for token_template_label, weight in weighted_entries:
                    if 'From training reaction' in token_template_label:
                        token_template_label = token_template_label.split()[-1]
                    rule_entry, training_entry = self.retrieve_original_entry(token_template_label)
                    if training_entry:
                        if (rule_entry, training_entry) in training:
                            training[(rule_entry, training_entry)] += weight
                        else:
                            training[(rule_entry, training_entry)] = weight
                    else:
                        if rule_entry in rules:
                            rules[rule_entry] += weight
                        else:
                            rules[rule_entry] = weight
                # Each entry should now only appear once
                training = [(k[0], k[1], v) for k, v in training.items()]
                rules = list(rules.items())

            return rules, training

    def extract_source_from_comments(self, reaction):
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
        lines = reaction.kinetics.comment.split('\n')

        exact_rule = False
        template = None
        rules = None
        training_entries = None
        degeneracy = 1

        training_reaction_pattern = r'Matched reaction\s*(\d+).*in.*training'
        degeneracy_pattern = r'Multiplied by reaction path degeneracy\s*(\d+)'

        for line in lines:
            training_matches = re.search(training_reaction_pattern, line)
            degeneracy_matches = re.search(degeneracy_pattern, line)

            if training_matches is not None:
                # Source of the kinetics is from training reaction
                training_reaction_index = int(training_matches.group(1))
                depository = self.get_training_depository()
                training_entry = depository.entries[training_reaction_index]
                # Perform sanity check that the training reaction's label matches that of the comments
                if training_entry.label not in line:
                    raise AssertionError(f'Reaction {reaction} uses kinetics from training reaction {training_reaction_index} '
                                         f'but does not match the training reaction {training_reaction_index} from the '
                                         f'{self.label} family.')

                # Sometimes the matched kinetics could be in the reverse direction.....
                if reaction.is_isomorphic(training_entry.item, either_direction=False, save_order=self.save_order):
                    reverse = False
                else:
                    reverse = True
                return True, [self.label, training_entry, reverse]

            if 'Exact match found for rate rule' in line:
                exact_rule = True
            if degeneracy_matches is not None:
                degeneracy = float(degeneracy_matches.group(1))

        # Extract the rate rule information
        full_comment_string = reaction.kinetics.comment.replace('\n', ' ')
        autogen_node_search_pattern = r'Estimated from node (.*)'
        # The rate rule string is right after the phrase 'for rate rule'
        template_pattern = r"for rate rule \[(.*)\]"  # only hit outermost brackets
        autogen_node_matches = re.search(autogen_node_search_pattern, full_comment_string)
        template_matches = re.search(template_pattern, full_comment_string)
        if autogen_node_matches is not None:  # autogenerated trees
            template_str = autogen_node_matches.group(1).split('Multiplied by reaction path degeneracy')[0].strip()
            template_str = template_str.split('in family')[0].strip()
            tokens = template_str.split()
            if len(tokens) == 2:  # The node was probably split because wordwrap was turned off
                assert len(template_str) > 115, 'The node name is too short to have been broken up by the chemkin writer'
                template_str = ''.join(tokens)
            elif len(tokens) > 2:  # warn the user the node is probably wrong
                raise ValueError(f'The node name {template_str} has multiple spaces and cannot be parsed for reaction {reaction}.')
            template = self.retrieve_template([template_str])
        elif template_matches is not None:  # hand-built trees
            template_label = template_matches.group(1)
            template = self.retrieve_template(template_label.split(';'))
        else:
            raise ValueError(f'Could not find rate rule in comments for reaction {reaction}.')
        rules, training_entries = self.get_sources_for_template(template)
        source_dict = {'template': template, 'degeneracy': degeneracy, 'exact': exact_rule,
                       'rules': rules, 'training': training_entries}

        # Source of the kinetics is from rate rules
        return False, [self.label, source_dict]

    def get_backbone_roots(self):
        """
        Returns: the top level backbone node in a unimolecular family.
        """

        backbone_roots = [entry for entry in self.groups.top if entry in self.forward_template.reactants]
        return backbone_roots

    def get_end_roots(self):
        """
        Returns: A list of top level end nodes in a unimolecular family
        """

        end_roots = [entry for entry in self.groups.top if entry not in self.forward_template.reactants]
        return end_roots

    def get_top_level_groups(self, root):
        """
        Returns a list of group nodes that are the highest in the tree starting at node "root".
        If "root" is a group node, then it will return a single-element list with "root".
        Otherwise, for every child of root, we descend until we find no nodes with logic
        nodes. We then return a list of all group nodes found along the way.
        """

        group_list = [root]
        all_groups = False

        while not all_groups:
            new_group_list = []
            for entry in group_list:
                if isinstance(entry.item, Group):
                    new_group_list.append(entry)
                else:
                    new_group_list.extend(entry.children)
            group_list = new_group_list
            all_groups = all([isinstance(entry.item, Group) for entry in group_list])

        return group_list


def information_gain(ks1, ks2):
    """
    calculates the information gain as the sum of the products of the standard deviations at each
    node and the number of reactions at that node
    """
    return len(ks1) * np.std(ks1) + len(ks2) * np.std(ks2)


def get_objective_function(kinetics1, kinetics2, obj=information_gain, T=1000.0):
    """
    Returns the value of four potential objective functions to minimize
    Uncertainty = N1*std(Ln(k))_1 + N1*std(Ln(k))_1
    Mean difference: -abs(mean(Ln(k))_1-mean(Ln(k))_2)
    Error using mean: Err_1 + Err_2
    Split: abs(N1-N2)
    """
    if not isinstance(kinetics1[0], Marcus):
        ks1 = np.array([np.log(k.get_rate_coefficient(T)) for k in kinetics1])
        ks2 = np.array([np.log(k.get_rate_coefficient(T)) for k in kinetics2])
    else:
        ks1 = np.array([k.get_lmbd_i(T) for k in kinetics1])
        ks2 = np.array([k.get_lmbd_i(T) for k in kinetics2])
    N1 = len(ks1)

    return obj(ks1, ks2), N1 == 0


def _make_rule(rr):
    """
    Function for parallelization of rule and uncertainty calculation

    Input: rr - tuple of (recipe, rxns, Tref, fmax, label, ranks)
           rxns and ranks are lists of equal length.
    Output: kinetics object, with uncertainty and comment attached.
    If Blowers-Masel fitting is successful it will be ArrheniusBM or ArrheniusChargeTransferBM,
    else Arrhenius, SurfaceChargeTransfer, or ArrheniusChargeTransfer.

    Errors in Ln(k) at each reaction are treated as samples from a weighted normal distribution
    weights are inverse variance weights based on estimates of the error in Ln(k) for each individual reaction
    """
    recipe, rxns, Tref, fmax, label, ranks = rr
    for i, rxn in enumerate(rxns):
        rxn.rank = ranks[i]
    rxns = np.array(rxns)
    rs = np.array([r for r in rxns if type(r.kinetics) is not KineticsModel]) # KineticsModel is the base class with no data.
    n = len(rs)
    if n == 0:
        return None

    if isinstance(rs[0].kinetics, Marcus):
        kin = average_kinetics([r.kinetics for r in rs])
        return kin

    data_mean = np.mean(np.log([r.kinetics.get_rate_coefficient(Tref) for r in rs]))

    if isinstance(rs[0].kinetics, Arrhenius):
        arr = ArrheniusBM
    else:
        arr = ArrheniusChargeTransferBM
    if n > 1:
        kin = arr().fit_to_reactions(rs, recipe=recipe)
    if n == 1 or kin.E0.value_si < 0.0 or abs(kin.n.value_si) > 5.0:
        # still run it through the averaging function when n=1 to standardize the units and run checks
        kin = average_kinetics([r.kinetics for r in rs])
        if n == 1:
            kin.uncertainty = RateUncertainty(mu=0.0, var=(np.log(fmax) / 2.0) ** 2, N=1, Tref=Tref, data_mean=data_mean, correlation=label)
            kin.comment = f"Only one reaction rate: {rs[0]!s}"
        else:
            if kin.E0.value_si < 0.0:
                reason = "E0<0"
            elif abs(kin.n.value_si) > 5.0:
                reason = "n>5"
            else:
                reason = "?"
            kin.comment = f"Blowers-Masel fit was bad ({reason}) so instead averaged from {n} reactions."
            dlnks = np.array([
                np.log(
                        average_kinetics([r.kinetics for r in rs[list(set(range(len(rs))) - {i})]]).get_rate_coefficient(T=Tref) / rxn.get_rate_coefficient(T=Tref)
                    ) for i, rxn in enumerate(rs)
                ])  # 1) fit to set of reactions without the current reaction (k)  2) compute log(kfit/kactual) at Tref
            varis = (np.array([rank_accuracy_map[rxn.rank].value_si for rxn in rs]) / (2.0 * 8.314 * Tref)) ** 2
            # weighted average calculations
            ws = 1.0 / varis
            V1 = ws.sum()
            V2 = (ws ** 2).sum()
            mu = np.dot(ws, dlnks) / V1
            s = np.sqrt(np.dot(ws, (dlnks - mu) ** 2) / (V1 - V2 / V1))
            kin.uncertainty = RateUncertainty(mu=mu, var=s ** 2, N=n, Tref=Tref, data_mean=data_mean, correlation=label)
    else: # Blowers-Masel fit was good
        if isinstance(rs[0].kinetics, Arrhenius):
            dlnks = np.array([
                np.log(
                    arr().fit_to_reactions(rs[list(set(range(len(rs))) - {i})], recipe=recipe)
                    .to_arrhenius(rxn.get_enthalpy_of_reaction(298.))
                    .get_rate_coefficient(T=Tref) / rxn.get_rate_coefficient(T=Tref)
                ) for i, rxn in enumerate(rs)
            ])  # 1) fit to set of reactions without the current reaction (k)  2) compute log(kfit/kactual) at Tref
        else: # SurfaceChargeTransfer or ArrheniusChargeTransfer
            dlnks = np.array([
                np.log(
                    arr().fit_to_reactions(rs[list(set(range(len(rs))) - {i})], recipe=recipe)
                    .to_arrhenius_charge_transfer(rxn.get_enthalpy_of_reaction(298.))
                    .get_rate_coefficient(T=Tref) / rxn.get_rate_coefficient(T=Tref)
                ) for i, rxn in enumerate(rs)
            ])  # 1) fit to set of reactions without the current reaction (k)  2) compute log(kfit/kactual) at Tref
        varis = (np.array([rank_accuracy_map[rxn.rank].value_si for rxn in rs]) / (2.0 * 8.314 * Tref)) ** 2
        # weighted average calculations
        ws = 1.0 / varis
        V1 = ws.sum()
        V2 = (ws ** 2).sum()
        mu = np.dot(ws, dlnks) / V1
        s = np.sqrt(np.dot(ws, (dlnks - mu) ** 2) / (V1 - V2 / V1))
        kin.uncertainty = RateUncertainty(mu=mu, var=s ** 2, N=n, Tref=Tref, data_mean=data_mean, correlation=label)

    #site solute parameters
    site_datas = [get_site_solute_data(rxn) for rxn in rxns]
    site_datas = [sdata for sdata in site_datas if sdata is not None]
    if len(site_datas) > 0:
        site_data = SoluteTSData()
        for sdata in site_datas:
            site_data += sdata
        site_data = site_data * (1.0/len(site_datas))
        kin.solute = site_data
    return kin

def _spawn_tree_process(family, template_rxn_map, obj, T, nprocs, depth, min_splitable_entry_num, min_rxns_to_spawn, extension_iter_max, extension_iter_item_cap):
    parent_conn, child_conn = mp.Pipe()
    name = list(template_rxn_map.keys())[0]
    p = mp.Process(target=_child_make_tree_nodes,
                   args=(family, child_conn, template_rxn_map, obj, T, nprocs,
                         depth, name, min_splitable_entry_num, min_rxns_to_spawn, extension_iter_max, extension_iter_item_cap))
    p.start()
    return parent_conn, p, name


def _child_make_tree_nodes(family, child_conn, template_rxn_map, obj, T, nprocs, depth, name, min_splitable_entry_num,
                           min_rxns_to_spawn, extension_iter_max, extension_iter_item_cap):
    del_labels = []
    root_label = list(template_rxn_map.keys())[0]
    for label in family.groups.entries.keys():
        if label != root_label:
            del_labels.append(label)
    for label in del_labels:
        del family.groups.entries[label]

    family.groups.entries[root_label].parent = None

    family.make_tree_nodes(template_rxn_map=template_rxn_map, obj=obj, T=T, nprocs=nprocs, depth=depth + 1,
                           min_splitable_entry_num=min_splitable_entry_num, min_rxns_to_spawn=min_rxns_to_spawn,
                           extension_iter_max=extension_iter_max, extension_iter_item_cap=extension_iter_item_cap)

    child_conn.send(list(family.groups.entries.values()))

def average_kinetics(kinetics_list):
    """
    Based on averaging log k.
    Hence we average n, Ea, arithmetically, but we
    average log A (geometric average)
    """
    if type(kinetics_list[0]) not in [Arrhenius,SurfaceChargeTransfer,ArrheniusChargeTransfer,Marcus]:
        raise Exception('Invalid kinetics type {0!r} for {1!r}.'.format(type(kinetics), self))
    
    Aunits = kinetics_list[0].A.units
    if Aunits in {'cm^3/(mol*s)', 'cm^3/(molecule*s)', 'm^3/(molecule*s)'}:
        Aunits = 'm^3/(mol*s)'
    elif Aunits in {'cm^6/(mol^2*s)', 'cm^6/(molecule^2*s)', 'm^6/(molecule^2*s)'}:
        Aunits = 'm^6/(mol^2*s)'
    elif Aunits in {'s^-1', 'm^3/(mol*s)', 'm^6/(mol^2*s)'}:
        # they were already in SI
        pass
    elif Aunits in {'m^2/(mol*s)', 'cm^2/(mol*s)', 'm^2/(molecule*s)', 'cm^2/(molecule*s)'}:
        # surface: bimolecular (Langmuir-Hinshelwood)
        Aunits = 'm^2/(mol*s)'
    elif Aunits in {'m^5/(mol^2*s)', 'cm^5/(mol^2*s)', 'm^5/(molecule^2*s)', 'cm^5/(molecule^2*s)'}:
        # surface: dissociative adsorption
        Aunits = 'm^5/(mol^2*s)'
    elif Aunits == '':
        # surface: sticking coefficient
        pass
    else:
        raise Exception('Invalid units {0} for averaging kinetics.'.format(Aunits))
    
    logA = 0.0
    n = 0.0
    Ea = 0.0
    alpha = 0.5
    lmbd_i_coefs = np.zeros(4)
    beta = 0.0
    wr = 0.0
    wp = 0.0
    electrons = None
    if isinstance(kinetics_list[0], SurfaceChargeTransfer) or isinstance(kinetics_list[0], ArrheniusChargeTransfer):
        if electrons is None:
            electrons = kinetics_list[0].electrons.value_si
        assert all(np.abs(k.V0.value_si) < 0.0001 for k in kinetics_list), [k.V0.value_si for k in kinetics_list]
        assert all(np.abs(k.alpha.value_si - 0.5) < 0.001 for k in kinetics_list), [k.alpha for k in kinetics_list]
    V0 = 0.0
    count = 0
    for kinetics in kinetics_list:
        count += 1
        logA += np.log10(kinetics.A.value_si)
        n += kinetics.n.value_si
        if hasattr(kinetics,"Ea"):
            Ea += kinetics.Ea.value_si
        if hasattr(kinetics,"lmbd_i_coefs"):
            lmbd_i_coefs += kinetics.lmbd_i_coefs.value_si
            beta += kinetics.beta.value_si
            wr += kinetics.wr.value_si
            wp += kinetics.wp.value_si

    logA /= count
    n /= count
    Ea /= count
    lmbd_i_coefs /= count
    beta /= count 
    wr /= count 
    wp /= count

    if isinstance(kinetics, Marcus):
        averaged_kinetics = Marcus(
            A=(10 ** logA, Aunits),
            n=n,
            lmbd_i_coefs=lmbd_i_coefs,
            beta=(beta,"1/m"),
            wr=(wr * 0.001, "kJ/mol"),
            wp=(wp * 0.001, "kJ/mol"),
            comment=f"Averaged from {len(kinetics_list)} rate expressions.",
            )
    elif isinstance(kinetics, SurfaceChargeTransfer):
        averaged_kinetics = SurfaceChargeTransfer(
            A=(10 ** logA, Aunits),
            n=n,
            electrons=electrons,
            alpha=alpha,
            V0=(V0,'V'),
            Ea=(Ea * 0.001, "kJ/mol"),
            comment=f"Averaged from {len(kinetics_list)} rate expressions.",
            )
    elif isinstance(kinetics, ArrheniusChargeTransfer):
        averaged_kinetics = ArrheniusChargeTransfer(
            A=(10 ** logA, Aunits),
            n=n,
            electrons=electrons,
            alpha=alpha,
            V0=(V0,'V'),
            Ea=(Ea * 0.001, "kJ/mol"),
            comment=f"Averaged from {len(kinetics_list)} rate expressions.",
            )
    else:
        averaged_kinetics = Arrhenius(
            A=(10 ** logA, Aunits),
            n=n,
            Ea=(Ea * 0.001, "kJ/mol"),
            comment=f"Averaged from {len(kinetics_list)} rate expressions.",
        )
    return averaged_kinetics

def get_site_solute_data(rxn):
    """
    apply kinetic solvent correction in this case the parameters are dGTSsite instead of GTS
    """
    from rmgpy.data.rmg import get_db
    solvation_database = get_db('solvation')
    ts_data = rxn.kinetics.solute
    if ts_data:
        site_data = to_soluteTSdata(ts_data,reactants=rxn.reactants)

        #compute x from gas phase
        GR = 0.0
        GP = 0.0

        for reactant in rxn.reactants:
            try:
                GR += reactant.thermo.get_free_energy(298.0)
            except Exception:
                logging.error("Problem with reactant {!r} in reaction {!s}".format(reactant, rxn))
                raise
        for product in rxn.products:
            try:
                GP += product.thermo.get_free_energy(298.0)
            except Exception:
                logging.error("Problem with product {!r} in reaction {!s}".format(reactant, rxn))
                raise

        dGrxn = GP-GR
        if dGrxn > 0:
            x = 1.0
        else:
            x = 0.0

        for spc in rxn.reactants:
            spc_solute_data = to_soluteTSdata(solvation_database.get_solute_data(spc.copy(deep=True)))
            site_data -= spc_solute_data*(1.0-x)

        for spc in rxn.products:
            spc_solute_data = to_soluteTSdata(solvation_database.get_solute_data(spc.copy(deep=True)))
            site_data -= spc_solute_data*x

        return site_data
    else:
        return None
