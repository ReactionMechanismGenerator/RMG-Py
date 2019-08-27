#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
This module provides classes and methods for working with molecular substructure
groups. These enable molecules to be searched for common motifs (e.g.
reaction sites).
"""

import itertools
from copy import deepcopy, copy

import cython

import rmgpy.molecule.element as elements
import rmgpy.molecule.molecule as mol
from rmgpy.exceptions import ActionError, ImplicitBenzeneError, UnexpectedChargeError
from rmgpy.molecule.atomtype import atomTypes, allElements, nonSpecifics, getFeatures, AtomType
from rmgpy.molecule.element import PeriodicSystem
from rmgpy.molecule.graph import Vertex, Edge, Graph


################################################################################

class GroupAtom(Vertex):
    """
    An atom group. This class is based on the :class:`Atom` class, except that
    it uses :ref:`atom types <atom-types>` instead of elements, and all
    attributes are lists rather than individual values. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atomType`          ``list``            The allowed atom types (as :class:`AtomType` objects)
    `radicalElectrons`  ``list``            The allowed numbers of radical electrons (as short integers)
    `charge`            ``list``            The allowed formal charges (as short integers)
    `label`             ``str``             A string label that can be used to tag individual atoms
    `lonePairs`         ``list``            The number of lone electron pairs
    'charge'            ''list''            The partial charge of the atom
    `props`             ``dict``            Dictionary for storing additional atom properties
    `reg_dim_atm`       ``list``            List of atom types that are free dimensions in tree optimization
    `reg_dim_u`         ``list``            List of unpaired electron numbers that are free dimensions in tree optimization
    `reg_dim_r`         ``list``            List of inRing values that are free dimensions in tree optimization
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. an atom will match the
    group if it matches *any* item in the list. However, the
    `radicalElectrons`, and `charge` attributes are linked
    such that an atom must match values from the same index in each of these in
    order to match.
    """

    def __init__(self, atomType=None, radicalElectrons=None, charge=None, label='', lonePairs=None, props=None):
        Vertex.__init__(self)
        self.atomType = atomType or []
        for index in range(len(self.atomType)):
            if isinstance(self.atomType[index], str):
                self.atomType[index] = atomTypes[self.atomType[index]]
        self.radicalElectrons = radicalElectrons or []
        self.charge = charge or []
        self.label = label
        self.lonePairs = lonePairs or []

        self.props = props or {}

        self.reg_dim_atm = [[], []]
        self.reg_dim_u = [[], []]
        self.reg_dim_r = [[], []]

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'edges': self.edges,
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sortingLabel': self.sortingLabel,
        }
        atomType = self.atomType
        if atomType is not None:
            atomType = [a.label for a in atomType]
        return (GroupAtom, (atomType, self.radicalElectrons, self.charge, self.label, self.lonePairs), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sortingLabel = d['sortingLabel']

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '[{0} {1}]'.format(self.label, ','.join([repr(a.label) for a in self.atomType]))

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<GroupAtom {0!s}>".format(self)

    @property
    def bonds(self):
        return self.edges

    def copy(self):
        """
        Return a deep copy of the :class:`GroupAtom` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupAtom(
            self.atomType[:],
            self.radicalElectrons[:],
            self.charge[:],
            self.label,
            self.lonePairs[:],
            deepcopy(self.props),
        )

    def __changeBond(self, order):
        """
        Update the atom group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        atomType = []
        for atom in self.atomType:
            if order == 1:
                atomType.extend(atom.incrementBond)
            elif order == -1:
                atomType.extend(atom.decrementBond)
            else:
                raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: '
                                  'Invalid order "{0}".'.format(order))
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __formBond(self, order):
        """
        Update the atom group as a result of applying a FORM_BOND action,
        where `order` specifies the order of the forming bond, and should be
        1 (since we only allow forming of single bonds).
        """
        if order == 0:
            # no change to atom types!
            return
        if order != 1:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.formBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __breakBond(self, order):
        """
        Update the atom group as a result of applying a BREAK_BOND action,
        where `order` specifies the order of the breaking bond, and should be
        1 (since we only allow breaking of single bonds).
        """
        if order == 0:
            # no change to atom types!
            return
        if order != 1:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.breakBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __gainRadical(self, radical):
        """
        Update the atom group as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.

        The 'radicalElectron' attribute can be an empty list if we use the wildcard
        argument ux in the group definition. In this case, we will have this
        function set the atom's 'radicalElectron' to a list allowing 1, 2, 3,
        or 4 radical electrons.
        """
        radicalElectrons = []
        if any([len(atomType.incrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_RADICAL action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))
        if not self.radicalElectrons:
            radicalElectrons = [1, 2, 3, 4]
        else:
            for electron in self.radicalElectrons:
                radicalElectrons.append(electron + radical)
        # Set the new radical electron counts
        self.radicalElectrons = radicalElectrons

    def __loseRadical(self, radical):
        """
        Update the atom group as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.

        The 'radicalElectron' attribute can be an empty list if we use the wildcard
        argument ux in the group definition. In this case, we will have this
        function set the atom's 'radicalElectron' to a list allowing 0, 1, 2,
        or 3 radical electrons.
        """
        radicalElectrons = []
        if any([len(atomType.decrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))

        if not self.radicalElectrons:
            radicalElectrons = [0, 1, 2, 3]
        else:
            for electron in self.radicalElectrons:
                electron = electron - radical
                if electron < 0:
                    raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: '
                                      'Invalid radical electron set "{0}".'.format(self.radicalElectrons))
                radicalElectrons.append(electron)

        # Set the new radical electron counts
        self.radicalElectrons = radicalElectrons

    def __gainPair(self, pair):
        """
        Update the atom group as a result of applying a GAIN_PAIR action,
        where `pair` specifies the number of lone electron pairs to add.
        """
        lonePairs = []
        atomType = []

        for atom in self.atomType:
            atomType.extend(atom.incrementLonePair)
        if any([len(atom.incrementLonePair) == 0 for atom in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_PAIR action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))

        # Add a lone pair to a group atom with none
        if not self.lonePairs:
            self.lonePairs = [1, 2, 3, 4]  # set to a wildcard of any number greater than 0
        # Add a lone pair to a group atom that already has at least one lone pair
        else:
            for x in self.lonePairs:
                lonePairs.append(x + pair)
            # Set the new lone electron pair count
            self.lonePairs = lonePairs

        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __losePair(self, pair):
        """
        Update the atom group as a result of applying a LOSE_PAIR action,
        where `pair` specifies the number of lone electron pairs to remove.
        """
        lonePairs = []
        atomType = []

        for atom in self.atomType:
            atomType.extend(atom.decrementLonePair)
        if any([len(atom.decrementLonePair) == 0 for atom in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomType))

        if not self.lonePairs:
            self.lonePairs = [0, 1, 2, 3]  # set to a wildcard of any number fewer than 4
        else:
            for x in self.lonePairs:
                if x - pair < 0:
                    raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: '
                                      'Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
                lonePairs.append(x - pair)
            # Set the new lone electron pair count
            self.lonePairs = lonePairs

        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def applyAction(self, action):
        """
        Update the atom group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        act = action[0].upper()
        if act == 'CHANGE_BOND':
            self.__changeBond(action[2])
        elif act == 'FORM_BOND':
            self.__formBond(action[2])
        elif act == 'BREAK_BOND':
            self.__breakBond(action[2])
        elif act == 'GAIN_RADICAL':
            self.__gainRadical(action[2])
        elif act == 'LOSE_RADICAL':
            self.__loseRadical(action[2])
        elif action[0].upper() == 'GAIN_PAIR':
            self.__gainPair(action[2])
        elif action[0].upper() == 'LOSE_PAIR':
            self.__losePair(action[2])
        else:
            raise ActionError('Unable to update GroupAtom: Invalid action {0}".'.format(action))

    def equivalent(self, other, strict=True):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Atom` or an :class:`GroupAtom`
        object. When comparing two :class:`GroupAtom` objects, this function
        respects wildcards, e.g. ``R!H`` is equivalent to ``C``.
        
        """
        cython.declare(group=GroupAtom)
        if not strict:
            raise NotImplementedError('There is currently no implementation of the strict argument for Group objects.')
        if not isinstance(other, GroupAtom):
            # Let the equivalent method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        group = other

        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short,
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType:
            for atomType2 in group.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        for atomType1 in group.atomType:
            for atomType2 in self.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        for radical1 in self.radicalElectrons:
            if group.radicalElectrons:  # Only check if the list is non-empty.  An empty list indicates a wildcard.
                for radical2 in group.radicalElectrons:
                    if radical1 == radical2: break
                else:
                    return False
        for radical1 in group.radicalElectrons:
            if self.radicalElectrons:
                for radical2 in self.radicalElectrons:
                    if radical1 == radical2: break
                else:
                    return False
        for lp1 in self.lonePairs:
            if group.lonePairs:
                for lp2 in group.lonePairs:
                    if lp1 == lp2: break
                else:
                    return False
        # Each charge in self must have an equivalent in other (and vice versa)
        for charge1 in self.charge:
            if group.charge:
                for charge2 in group.charge:
                    if charge1 == charge2: break
                else:
                    return False
        for charge1 in group.charge:
            if self.charge:
                for charge2 in self.charge:
                    if charge1 == charge2: break
                else:
                    return False
        # Other properties must have an equivalent in other (and vice versa)
        # Absence of the 'inRing' prop indicates a wildcard
        if 'inRing' in self.props and 'inRing' in group.props:
            if self.props['inRing'] != group.props['inRing']:
                return False
        # Otherwise the two atom groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `self` is the same as `other` or is a more
        specific case of `other`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive. 
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)
        group = other

        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short,
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType:  # all these must match
            for atomType2 in group.atomType:  # can match any of these
                if atomType1.isSpecificCaseOf(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        if self.radicalElectrons:
            for radical1 in self.radicalElectrons:
                if group.radicalElectrons:
                    for radical2 in group.radicalElectrons:
                        if radical1 == radical2: break
                    else:
                        return False
        else:
            if group.radicalElectrons: return False
        if self.lonePairs:
            for lp1 in self.lonePairs:
                if group.lonePairs:
                    for lp2 in group.lonePairs:
                        if lp1 == lp2: break
                    else:
                        return False
        else:
            if group.lonePairs: return False
        # Each charge in self must have an equivalent in other
        if self.charge:
            for charge1 in self.charge:
                if group.charge:
                    for charge2 in group.charge:
                        if charge1 == charge2: break
                    else:
                        return False
        else:
            if group.charge: return False
        # Other properties must have an equivalent in other
        # Absence of the 'inRing' prop indicates a wildcard
        if 'inRing' in self.props and 'inRing' in group.props:
            if self.props['inRing'] != group.props['inRing']:
                return False
        elif 'inRing' not in self.props and 'inRing' in group.props:
            return False
        # Otherwise self is in fact a specific case of other
        return True

    def isSurfaceSite(self):
        """
        Return ``True`` if the atom represents a surface site or ``False`` if not.
        """
        site_type = atomTypes['X']
        return all([s.isSpecificCaseOf(site_type) for s in self.atomType])

    def isOxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if not.
        """
        all_oxygens = [atomTypes['O']] + atomTypes['O'].specific
        check_list = [x in all_oxygens for x in self.atomType]
        return all(check_list)

    def isSulfur(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if not.
        """
        all_sulfur = [atomTypes['S']] + atomTypes['S'].specific
        check_list = [x in all_sulfur for x in self.atomType]
        return all(check_list)

    def isNitrogen(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if not.
        """
        all_nitrogen = [atomTypes['N']] + atomTypes['N'].specific
        check_list = [x in all_nitrogen for x in self.atomType]
        return all(check_list)

    def isCarbon(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if not.
        """
        all_carbon = [atomTypes['C']] + atomTypes['C'].specific
        check_list = [x in all_carbon for x in self.atomType]
        return all(check_list)

    def hasWildcards(self):
        """
        Return ``True`` if the atom has wildcards in any of the attributes:
        atomtype, radical electrons, lone pairs, charge, and bond order. Returns
        ''False'' if no attribute has wildcards.
        """
        if len(self.atomType) > 1:
            return True
        elif len(self.radicalElectrons) > 1 or len(self.radicalElectrons) == 0:
            return True
        elif len(self.lonePairs) > 1:
            return True
        for bond in self.bonds.values():
            if len(bond.order) > 1:
                return True
        return False

    def countBonds(self, wildcards=False):
        """
        Returns: list of the number of bonds currently on the :class:GroupAtom

        If the argument wildcards is turned off then any bonds with multiple
        options for bond orders will not be counted
        """
        # count up number of bonds
        single = r_double = o_double = s_double = triple = quadruple = benzene = 0  # note that 0 is immutable
        for atom2, bond12 in self.bonds.items():
            if not wildcards and len(bond12.order) > 1:
                continue
            # Count numbers of each higher-order bond type
            if bond12.isSingle(wildcards=True):
                single += 1
            if bond12.isDouble(wildcards=True):
                if atom2.isOxygen():
                    o_double += 1
                elif atom2.isSulfur():
                    s_double += 1
                else:
                    # r_double is for double bonds NOT to oxygen or Sulfur
                    r_double += 1
            if bond12.isTriple(wildcards=True): triple += 1
            if bond12.isQuadruple(wildcards=True): quadruple += 1
            if bond12.isBenzene(wildcards=True): benzene += 1

        all_double = r_double + o_double + s_double

        # Warning: some parts of code assume this matches precisely the list returned by getFeatures()
        return [single, all_double, r_double, o_double, s_double, triple, quadruple, benzene]

    def makeSampleAtom(self):
        """

        Returns: a class :Atom: object analagous to the GroupAtom

        This makes a sample, so it takes the first element when there are multiple options inside of
        self.atomtype, self.radicalElectrons, self.lonePairs, and self.charge

        """

        # Use the first atomtype to determine element, even if there is more than one atomtype
        atomtype = self.atomType[0]
        element = None

        default_lone_pairs = {'H': 0,
                              'D': 0,
                              'T': 0,
                              'He': 1,
                              'C': 0,
                              'O': 2,
                              'N': 1,
                              'Si': 0,
                              'S': 2,
                              'Ne': 4,
                              'Cl': 3,
                              'F': 3,
                              'I': 3,
                              'Ar': 4,
                              'X': 0,
                              }

        for element_label in allElements:
            if atomtype is atomTypes[element_label] or atomtype in atomTypes[element_label].specific:
                element = element_label
                break
        else:
            # For types that correspond to more than one type of element, pick the first that appears in specific
            for subtype in atomtype.specific:
                if subtype.label in allElements:
                    element = subtype.label
                    break

        # dummy default_atom to get default values
        default_atom = mol.Atom()

        # Three possible values for charge and lonePairs
        if self.charge:
            new_charge = self.charge[0]
        elif atomtype.charge:
            new_charge = atomtype.charge[0]
        else:
            new_charge = default_atom.charge

        if self.lonePairs:
            new_lone_pairs = self.lonePairs[0]
        elif atomtype.lonePairs:
            new_lone_pairs = atomtype.lonePairs[0]
        else:
            new_lone_pairs = default_atom.lonePairs

        new_atom = mol.Atom(
            element=element,
            radicalElectrons=self.radicalElectrons[0] if self.radicalElectrons else default_atom.radicalElectrons,
            charge=new_charge,
            lonePairs=new_lone_pairs,
            label=self.label if self.label else default_atom.label
        )

        # For some reason the default when no lone pairs is set to -100,
        # Based on git history, it is probably because RDKit requires a number instead of None
        if new_atom.lonePairs == -100:
            new_atom.lonePairs = default_lone_pairs[new_atom.symbol]

        return new_atom


################################################################################

class GroupBond(Edge):
    """
    A bond group. This class is based on the :class:`Bond` class, except that
    all attributes are lists rather than individual values. The allowed bond
    types are given :ref:`here <bond-types>`. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``list``            The allowed bond orders (as character strings)
    `reg_dim`           ``Boolean``         Indicates if this is a regularization dimension during tree generation
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. a bond will match the
    group if it matches *any* item in the list.
    """

    def __init__(self, atom1, atom2, order=None):
        Edge.__init__(self, atom1, atom2)
        if order is not None and all([isinstance(oneOrder, str) for oneOrder in order]):
            self.setOrderStr(order)
        elif order is not None and any([isinstance(oneOrder, str) for oneOrder in order]):
            raise ActionError('order list given {} does not consist of only strings or only numbers'.format(order))
        else:
            self.order = order or []

        self.reg_dim = [[], []]

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return str(self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<GroupBond {0!r}>".format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (GroupBond, (self.vertex1, self.vertex2, self.order))

    def copy(self):
        """
        Return a deep copy of the :class:`GroupBond` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupBond(self.vertex1, self.vertex2, self.order[:])

    def getOrderStr(self):
        """
        returns a list of strings representing the bond order
        """
        values = []
        for value in self.order:
            if value == 1:
                values.append('S')
            elif value == 2:
                values.append('D')
            elif value == 3:
                values.append('T')
            elif value == 4:
                values.append('Q')
            elif value == 1.5:
                values.append('B')
            elif value == 0:
                values.append('vdW')
            elif value == 0.1:
                values.append('H')
            else:
                raise TypeError('Bond order number {} is not hardcoded as a string'.format(value))
        return values

    def setOrderStr(self, newOrder):
        """
        set the bond order using a valid bond-order character list
        """

        values = []
        for value in newOrder:
            if value == 'S':
                values.append(1)
            elif value == 'D':
                values.append(2)
            elif value == 'T':
                values.append(3)
            elif value == 'Q':
                values.append(4)
            elif value == 'vdW':
                values.append(0)
            elif value == 'B':
                values.append(1.5)
            elif value == 'H':
                values.append(0.1)
            else:
                # try to see if an float disguised as a string was input by mistake
                try:
                    values.append(float(value))
                except ValueError:
                    raise TypeError('Bond order {} is not hardcoded into this method'.format(value))
        self.order = values

    def getOrderNum(self):
        """
        returns the bond order as a list of numbers
        """
        return self.order

    def setOrderNum(self, newOrder):
        """
        change the bond order with a list of numbers
        """
        self.order = newOrder

    def isSingle(self, wildcards=False):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are single.
        
        NOTE: we can replace the absolute value relation with math.isclose when
        we swtich to python 3.5+
        """
        if wildcards:
            for order in self.order:
                if abs(order - 1) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 1) <= 1e-9 and len(self.order) == 1

    def isDouble(self, wildcards=False):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are double.
        """
        if wildcards:
            for order in self.order:
                if abs(order - 2) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 2) <= 1e-9 and len(self.order) == 1

    def isTriple(self, wildcards=False):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are triple.
        """
        if wildcards:
            for order in self.order:
                if abs(order - 3) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 3) <= 1e-9 and len(self.order) == 1

    def isQuadruple(self, wildcards=False):
        """
        Return ``True`` if the bond represents a quadruple bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are quadruple.
        """
        if wildcards:
            for order in self.order:
                if abs(order - 4) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 4) <= 1e-9 and len(self.order) == 1

    def isVanDerWaals(self, wildcards=False):
        """
        Return ``True`` if the bond represents a van der Waals bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are van der Waals.
        """
        if wildcards:
            for order in self.order:
                if abs(order[0]) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0]) <= 1e-9 and len(self.order) == 1

    def isBenzene(self, wildcards=False):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are benzene
        """
        if wildcards:
            for order in self.order:
                if abs(order - 1.5) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 1.5) <= 1e-9 and len(self.order) == 1

    def isHydrogenBond(self, wildcards=False):
        """
        Return ``True`` if the bond represents a hydrogen bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are hydrogen bonds.
        """
        if wildcards:
            for order in self.order:
                if abs(order) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0]) <= 1e-9 and len(self.order) == 1

    def __changeBond(self, order):
        """
        Update the bond group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order. `order` is normally 1 or -1, but can be any value
        """
        new_order = [value + order for value in self.order]
        if any([value < 0 or value > 4 for value in new_order]):
            raise ActionError('Unable to update Bond due to CHANGE_BOND action: '
                              'Invalid resulting order "{0}".'.format(new_order))
        # Change any modified benzene orders to the appropriate stable order
        new_order = set(new_order)
        if 0.5 in new_order:
            new_order.remove(0.5)
            new_order.add(1)
        if 2.5 in new_order:
            new_order.remove(2.5)
            new_order.add(2)
        # Allow formation of benzene bonds if a double bond can be formed
        if 2 in new_order:
            new_order.add(1.5)
        # Set the new bond orders
        self.order = list(new_order)

    def applyAction(self, action):
        """
        Update the bond group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self.__changeBond(action[2])
        else:
            raise ActionError('Unable to update GroupBond: Invalid action {0}".'.format(action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Bond` or an :class:`GroupBond`
        object.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the equivalent method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        gb = other

        cython.declare(order1=float, order2=float)
        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for order1 in self.order:
            for order2 in gb.order:
                if order1 == order2: break
            else:
                return False
        for order1 in gb.order:
            for order2 in self.order:
                if order1 == order2: break
            else:
                return False
        # Otherwise the two bond groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)
        gb = other

        cython.declare(order1=float, order2=float)
        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other
        for order1 in self.order:  # all these must match
            for order2 in gb.order:  # can match any of these
                if order1 == order2: break
            else:
                return False
        # Otherwise self is in fact a specific case of other
        return True

    def makeBond(self, molecule, atom1, atom2):
        """
        Creates a :class: Bond between atom1 and atom2 analogous to self

        The intended input arguments should be class :Atom: not class :GroupAtom:
        Args:
            atom1: First :class: Atom the bond connects
            atom2: Second :class: Atom the bond connects

        """
        new_bond = mol.Bond(atom1, atom2, order=self.order[0])
        molecule.addBond(new_bond)


class Group(Graph):
    """
    A representation of a molecular substructure group using a graph data
    type, extending the :class:`Graph` class. The attributes are:
    
    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atoms`             ``list``            Aliases for the `vertices` storing :class:`GroupAtom`
    `multiplicity`      ``list``            Range of multiplicities accepted for the group
    `props`             ``dict``            Dictionary of arbitrary properties/flags classifying state of Group object 
    =================== =================== ====================================

    Corresponding alias methods to Molecule have also been provided.
    """

    def __init__(self, atoms=None, props=None, multiplicity=None):
        Graph.__init__(self, atoms)
        self.props = props or {}
        self.multiplicity = multiplicity or []
        self.elementCount = {}
        self.radicalCount = -1
        self.update()

    def __deepcopy__(self, memo):
        return self.copy(deep=True)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Group, (self.vertices, self.props))

    def _repr_png_(self):
        """
        Return a png picture of the group, useful for ipython-qtconsole.
        """
        return self.draw('png')

    def draw(self, format):
        """
        Use pydot to draw a basic graph of the group.

        Use format to specify the desired output format, eg. 'png', 'svg', 'ps', 'pdf', 'plain', etc.
        """
        import pydot

        graph = pydot.Dot(graph_type='graph', dpi="52")
        for index, atom in enumerate(self.atoms):
            atom_type = '{0!s} '.format(atom.label if atom.label != '' else '')
            atom_type += ','.join([at.label for at in atom.atomType])
            atom_type = '"' + atom_type + '"'
            graph.add_node(pydot.Node(name=str(index + 1), label=atom_type, fontname="Helvetica", fontsize="16"))
        for atom1 in self.atoms:
            for atom2, bond in atom1.bonds.items():
                index1 = self.atoms.index(atom1)
                index2 = self.atoms.index(atom2)
                if index1 < index2:
                    bond_type = ','.join([order for order in bond.getOrderStr()])
                    bond_type = '"' + bond_type + '"'
                    graph.add_edge(pydot.Edge(src=str(index1 + 1), dst=str(index2 + 1),
                                              label=bond_type, fontname="Helvetica", fontsize="16"))

        img = graph.create(prog='neato', format=format)
        return img

    @property
    def atoms(self):
        """
        List of atoms contained in the current molecule.

        Renames the inherited vertices attribute of :class:`Graph`.
        """
        return self.vertices

    @atoms.setter
    def atoms(self, atoms):
        self.vertices = atoms

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.addVertex(atom)

    def addBond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        return self.addEdge(bond)

    def getBonds(self, atom):
        """
        Return a list of the bonds involving the specified `atom`.
        """
        return self.getEdges(atom)

    def getBond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.getEdge(atom1, atom2)

    def hasAtom(self, atom):
        """
        Returns ``True`` if `atom` is an atom in the graph, or ``False`` if
        not.
        """
        return self.hasVertex(atom)

    def hasBond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.hasEdge(atom1, atom2)

    def containsSurfaceSite(self):
        """
        Returns ``True`` iff the group contains an 'X' surface site.
        """
        cython.declare(atom=GroupAtom)
        for atom in self.atoms:
            if atom.isSurfaceSite():
                return True
        return False

    def isSurfaceSite(self):
        "Returns ``True`` iff the group is nothing but a surface site 'X'."
        return len(self.atoms) == 1 and self.atoms[0].isSurfaceSite()

    def removeAtom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        return self.removeVertex(atom)

    def removeBond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.removeEdge(bond)

    def removeVanDerWaalsBonds(self):
        """
        Remove all bonds that are definitely only van der Waals bonds.
        """
        cython.declare(atom=GroupAtom, bond=GroupBond)
        for bond in self.getAllEdges():
            if bond.isVanDerWaals(wildcards=False):
                self.removeBond(bond)

    def sortAtoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        return self.sortVertices()

    def sortByConnectivity(self, atomList):
        """
        Args:
            atomList: input list of atoms

        Returns: a sorted list of atoms where each atom is connected to a previous
        atom in the list if possible
        """
        # if no input given just return
        if not atomList: return atomList

        sorted_atom_list = []
        sorted_atom_list.append(atomList.pop(0))
        while atomList:
            for atom1 in sorted_atom_list:
                added = False
                for atom2, bond12 in atom1.bonds.items():
                    if bond12.isBenzene() and atom2 in atomList:
                        sorted_atom_list.append(atom2)
                        atomList.remove(atom2)
                        added = True
                if added: break
            else:
                sorted_atom_list.append(atomList.pop(0))

        return sorted_atom_list

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(Group)
        g = Graph.copy(self, deep)
        other = Group(g.vertices)
        return other

    def update(self):

        self.updateConnectivityValues()
        self.updateFingerprint()

    def update_charge(self):
        """
        Update the partial charge according to the valence electron, total bond order, lone pairs
        and radical electrons. This method is used for products of specific families with recipes that modify charges.
        """
        for atom in self.atoms:
            if (len(atom.charge) == 1) and (len(atom.lonePairs) == 1) and (len(atom.radicalElectrons) == 1):
                # if the charge of the group is not labeled, then no charge update will be
                # performed. If there multiple charges are assigned, no update either.
                # Besides, this groupatom should have enough information to be updated
                atom_type = atom.atomType[0]
                for element in allElements:
                    if atom_type is atomTypes[element] or atom_type in atomTypes[element].specific:
                        bond_order = 0
                        valence_electron = elements.PeriodicSystem.valence_electrons[element]
                        for _, bond in atom.bonds.items():
                            bond_order += bond.order[0]
                        lone_pairs = atom.lonePairs[0]
                        radical_electrons = atom.radicalElectrons[0]
                        atom.charge[0] = valence_electron - bond_order - 2 * lone_pairs - radical_electrons
                    else:
                        # if the group is not specified to specific element, charge will not be updated
                        pass

    def getNetCharge(self):
        """
        Iterate through the atoms in the group and calculate the net charge
        """
        return sum([atom.charge[0] for atom in self.vertices if atom.charge != []])

    def merge(self, other):
        """
        Merge two groups so as to store them in a single
        :class:`Group` object. The merged :class:`Group`
        object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Group(atoms=g.vertices)
        return molecule

    def split(self):
        """
        Convert a single :class:`Group` object containing two or more
        unconnected groups into separate class:`Group` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Group(atoms=g.vertices)
            molecules.append(molecule)
        return molecules

    def clearRegDims(self):
        """
        clear regularization dimensions
        """
        cython.declare(atm=GroupAtom)
        for atm in self.atoms:
            atm.reg_dim_atm = [[], []]
            atm.reg_dim_u = [[], []]
            atm.reg_dim_r = [[], []]
        for bd in self.getAllEdges():
            bd.reg_dim = [[], []]

    def getExtensions(self, R=None, basename='', atmInd=None, atmInd2=None, Nsplits=None):
        """
        generate all allowed group extensions and their complements
        note all atomtypes except for elements and R/R!H's must be removed
        """
        cython.declare(atoms=list, atm=GroupAtom, atm2=GroupAtom, bd=GroupBond, i=int, j=int,
                       extents=list, RnH=list, typ=list)

        extents = []

        if Nsplits is None:
            Nsplits = len(self.split())

        # generate appropriate R and R!H
        if R is None:
            R = elements.BDE_elements  # set of possible R elements/atoms
            R = [atomTypes[x] for x in R]

        Rbonds = [1, 2, 3, 1.5]
        Run = [0, 1, 2, 3]

        RnH = R[:]
        RnH.remove(atomTypes['H'])

        atoms = self.atoms
        if atmInd is None:
            for i, atm in enumerate(atoms):
                typ = atm.atomType
                if atm.reg_dim_atm[0] == []:
                    if len(typ) == 1:
                        if typ[0].label == 'R':
                            extents.extend(self.specifyAtomExtensions(i, basename, R))  # specify types of atoms
                        elif typ[0].label == 'R!H':
                            extents.extend(self.specifyAtomExtensions(i, basename, RnH))
                    else:
                        extents.extend(self.specifyAtomExtensions(i, basename, typ))
                else:
                    if len(typ) == 1:
                        if typ[0].label == 'R':
                            extents.extend(
                                self.specifyAtomExtensions(i, basename, atm.reg_dim_atm[0]))  # specify types of atoms
                        elif typ[0].label == 'R!H':
                            extents.extend(
                                self.specifyAtomExtensions(i, basename, list(set(atm.reg_dim_atm[0]) & set(R))))
                    else:
                        extents.extend(
                            self.specifyAtomExtensions(i, basename, list(set(typ) & set(atm.reg_dim_atm[0]))))
                if atm.reg_dim_u[0] == []:
                    if len(atm.radicalElectrons) != 1:
                        if len(atm.radicalElectrons) == 0:
                            extents.extend(self.specifyUnpairedExtensions(i, basename, Run))
                        else:
                            extents.extend(self.specifyUnpairedExtensions(i, basename, atm.radicalElectrons))
                else:
                    if len(atm.radicalElectrons) != 1 and len(atm.reg_dim_u[0]) != 1:
                        if len(atm.radicalElectrons) == 0:
                            extents.extend(self.specifyUnpairedExtensions(i, basename, atm.reg_dim_u[0]))
                        else:
                            extents.extend(self.specifyUnpairedExtensions(i, basename, list(
                                set(atm.radicalElectrons) & set(atm.reg_dim_u[0]))))
                if atm.reg_dim_r[0] == [] and 'inRing' not in atm.props:
                    extents.extend(self.specifyRingExtensions(i, basename))

                extents.extend(self.specifyExternalNewBondExtensions(i, basename, Rbonds))
                for j, atm2 in enumerate(atoms):
                    if j < i and not self.hasBond(atm, atm2):
                        extents.extend(self.specifyInternalNewBondExtensions(i, j, Nsplits, basename, Rbonds))
                    elif j < i:
                        bd = self.getBond(atm, atm2)
                        if len(bd.order) > 1 and bd.reg_dim[0] == []:
                            extents.extend(self.specifyBondExtensions(i, j, basename, bd.order))
                        elif len(bd.order) > 1 and len(bd.reg_dim[0]) > 1 and len(bd.reg_dim[0]) > len(bd.reg_dim[1]):
                            extents.extend(self.specifyBondExtensions(i, j, basename, bd.reg_dim[0]))

        elif atmInd is not None and atmInd2 is not None:  # if both atmInd and atmInd2 are defined only look at the bonds between them
            i = atmInd
            j = atmInd2
            atm = atoms[i]
            atm2 = atoms[j]
            if j < i and not self.hasBond(atm, atm2):
                extents.extend(self.specifyInternalNewBondExtensions(i, j, Nsplits, basename, Rbonds))
            if self.hasBond(atm, atm2):
                bd = self.getBond(atm, atm2)
                if len(bd.order) > 1 and bd.reg_dim[0] == []:
                    extents.extend(self.specifyBondExtensions(i, j, basename, bd.order))
                elif len(bd.order) > 1 and len(bd.reg_dim[0]) > 1 and len(bd.reg_dim[0]) > len(bd.reg_dim[1]):
                    extents.extend(self.specifyBondExtensions(i, j, basename, bd.reg_dim[0]))

        elif atmInd is not None:  # look at the atom at atmInd
            i = atmInd
            atm = atoms[i]
            typ = atm.atomType
            if atm.reg_dim_atm[0] == []:
                if len(typ) == 1:
                    if typ[0].label == 'R':
                        extents.extend(self.specifyAtomExtensions(i, basename, R))  # specify types of atoms
                    elif typ[0].label == 'R!H':
                        extents.extend(self.specifyAtomExtensions(i, basename, RnH))
                else:
                    extents.extend(self.specifyAtomExtensions(i, basename, typ))
            else:
                if len(typ) == 1:
                    if typ[0].label == 'R':
                        extents.extend(
                            self.specifyAtomExtensions(i, basename, atm.reg_dim_atm[0]))  # specify types of atoms
                    elif typ[0].label == 'R!H':
                        extents.extend(self.specifyAtomExtensions(i, basename, list(set(atm.reg_dim_atm[0]) & set(R))))
                else:
                    extents.extend(self.specifyAtomExtensions(i, basename, list(set(typ) & set(atm.reg_dim_atm[0]))))
            if atm.reg_dim_u == []:
                if len(atm.radicalElectrons) != 1:
                    if len(atm.radicalElectrons) == 0:
                        extents.extend(self.specifyUnpairedExtensions(i, basename, Run))
                    else:
                        extents.extend(self.specifyUnpairedExtensions(i, basename, atm.radicalElectrons))
            else:
                if len(atm.radicalElectrons) != 1 and len(atm.reg_dim_u[0]) != 1:
                    if len(atm.radicalElectrons) == 0:
                        extents.extend(self.specifyUnpairedExtensions(i, basename, atm.reg_dim_u[0]))
                    else:
                        extents.extend(self.specifyUnpairedExtensions(i, basename, list(
                            set(atm.radicalElectrons) & set(atm.reg_dim_u[0]))))
            if atm.reg_dim_r[0] == [] and 'inRing' not in atm.props:
                extents.extend(self.specifyRingExtensions(i, basename))

            extents.extend(self.specifyExternalNewBondExtensions(i, basename, Rbonds))
            for j, atm2 in enumerate(atoms):
                if j < i and not self.hasBond(atm, atm2):
                    extents.extend(self.specifyInternalNewBondExtensions(i, j, Nsplits, basename, Rbonds))
                elif j < i:
                    bd = self.getBond(atm, atm2)
                    if len(bd.order) > 1 and bd.reg_dim == []:
                        extents.extend(self.specifyBondExtensions(i, j, basename, bd.order))
                    elif len(bd.order) > 1 and len(bd.reg_dim[0]) > 1 and len(bd.reg_dim[0]) > len(bd.reg_dim[1]):
                        extents.extend(self.specifyBondExtensions(i, j, basename, bd.reg_dim[0]))

        else:
            raise ValueError('atmInd must be defined if atmInd2 is defined')

        return extents

    def specifyAtomExtensions(self, i, basename, R):
        """
        generates extensions for specification of the type of atom defined by a given atomtype
        or set of atomtypes
        """
        cython.declare(grps=list, labelList=list, Rset=set, item=AtomType, grp=Group, grpc=Group, k=AtomType, p=str)

        grps = []
        Rset = set(R)
        for item in R:
            grp = deepcopy(self)
            grpc = deepcopy(self)
            old_atom_type = grp.atoms[i].atomType
            grp.atoms[i].atomType = [item]
            grpc.atoms[i].atomType = list(Rset - {item})

            if len(old_atom_type) > 1:
                labelList = []
                old_atom_type_str = ''
                for k in old_atom_type:
                    labelList.append(k.label)
                for p in sorted(labelList):
                    old_atom_type_str += p
            else:
                old_atom_type_str = old_atom_type[0].label

            grps.append(
                (grp, grpc, basename + '_' + str(i + 1) + old_atom_type_str + '->' + item.label, 'atomExt', (i,)))

        return grps

    def specifyRingExtensions(self, i, basename):
        """
        generates extensions for specifying if an atom is in a ring
        """
        cython.declare(grps=list, label_list=list, grp=Group, grpc=Group, atom_type=list, atom_type_str=str, k=AtomType,
                       p=str)

        grps = []
        label_list = []

        grp = deepcopy(self)
        grpc = deepcopy(self)
        grp.atoms[i].props['inRing'] = True
        grpc.atoms[i].props['inRing'] = False

        atom_type = grp.atoms[i].atomType

        if len(atom_type) > 1:
            atom_type_str = ''
            for k in atom_type:
                label_list.append(k.label)
            for p in sorted(label_list):
                atom_type_str += p
        else:
            atom_type_str = atom_type[0].label

        grps.append((grp, grpc, basename + '_' + str(i + 1) + atom_type_str + '-inRing', 'ringExt', (i,)))

        return grps

    def specifyUnpairedExtensions(self, i, basename, Run):
        """
        generates extensions for specification of the number of electrons on a given atom
        """

        grps = []
        label_list = []

        Rset = set(Run)
        for item in Run:
            grp = deepcopy(self)
            grpc = deepcopy(self)
            grp.atoms[i].radicalElectrons = [item]
            grpc.atoms[i].radicalElectrons = list(Rset - {item})

            atom_type = grp.atoms[i].atomType

            if len(atom_type) > 1:
                atom_type_str = ''
                for k in atom_type:
                    label_list.append(k.label)
                for p in sorted(label_list):
                    atom_type_str += p
            else:
                atom_type_str = atom_type[0].label

            grps.append((grp, grpc, basename + '_' + str(i + 1) + atom_type_str + '-u' + str(item), 'elExt', (i,)))

        return grps

    def specifyInternalNewBondExtensions(self, i, j, Nsplits, basename, Rbonds):
        """
        generates extensions for creation of a bond (of undefined order)
        between two atoms indexed i,j that already exist in the group and are unbonded
        """
        cython.declare(newgrp=Group)

        label_list = []

        newgrp = deepcopy(self)
        newgrp.addBond(GroupBond(newgrp.atoms[i], newgrp.atoms[j], Rbonds))

        atom_type_i = newgrp.atoms[i].atomType
        atom_type_j = newgrp.atoms[j].atomType

        if len(atom_type_i) > 1:
            atom_type_i_str = ''
            for k in atom_type_i:
                label_list.append(k.label)
            for k in sorted(label_list):
                atom_type_i_str += k
        else:
            atom_type_i_str = atom_type_i[0].label
        if len(atom_type_j) > 1:
            atom_type_j_str = ''
            for k in atom_type_j:
                label_list.append(k.label)
            for p in sorted(label_list):
                atom_type_j_str += p
        else:
            atom_type_j_str = atom_type_j[0].label

        if len(newgrp.split()) < Nsplits:  # if this formed a bond between two seperate groups in the
            return []
        else:
            return [(newgrp, None,
                     basename + '_Int-' + str(i + 1) + atom_type_i_str + '-' + str(j + 1) + atom_type_j_str,
                     'intNewBondExt', (i, j))]

    def specifyExternalNewBondExtensions(self, i, basename, Rbonds):
        """
        generates extensions for the creation of a bond (of undefined order) between
        an atom and a new atom that is not H
        """
        cython.declare(ga=GroupAtom, newgrp=Group, j=int)

        label_list = []

        ga = GroupAtom([atomTypes['R!H']])
        newgrp = deepcopy(self)
        newgrp.addAtom(ga)
        j = newgrp.atoms.index(ga)
        newgrp.addBond(GroupBond(newgrp.atoms[i], newgrp.atoms[j], Rbonds))
        atom_type = newgrp.atoms[i].atomType
        if len(atom_type) > 1:
            atom_type_str = ''
            for k in atom_type:
                label_list.append(k.label)
            for p in sorted(label_list):
                atom_type_str += p
        else:
            atom_type_str = atom_type[0].label

        return [(newgrp, None, basename + '_Ext-' + str(i + 1) + atom_type_str + '-R', 'extNewBondExt',
                 (len(newgrp.atoms) - 1,))]

    def specifyBondExtensions(self, i, j, basename, Rbonds):
        """
        generates extensions for the specification of bond order for a given bond
        """
        cython.declare(grps=list, label_list=list, Rbset=set, bd=float, grp=Group, grpc=Group)
        grps = []
        label_list = []
        Rbset = set(Rbonds)
        bdict = {1: '-', 2: '=', 3: '#', 1.5: '-='}
        for bd in Rbonds:
            grp = deepcopy(self)
            grpc = deepcopy(self)
            grp.atoms[i].bonds[grp.atoms[j]].order = [bd]
            grp.atoms[j].bonds[grp.atoms[i]].order = [bd]
            grpc.atoms[i].bonds[grpc.atoms[j]].order = list(Rbset - {bd})
            grpc.atoms[j].bonds[grpc.atoms[i]].order = list(Rbset - {bd})

            atom_type_i = grp.atoms[i].atomType
            atom_type_j = grp.atoms[j].atomType

            if len(atom_type_i) > 1:
                atom_type_i_str = ''
                for k in atom_type_i:
                    label_list.append(k.label)
                for p in sorted(label_list):
                    atom_type_i_str += p
            else:
                atom_type_i_str = atom_type_i[0].label
            if len(atom_type_j) > 1:
                atom_type_j_str = ''
                for k in atom_type_j:
                    label_list.append(k.label)
                for p in sorted(label_list):
                    atom_type_j_str += p
            else:
                atom_type_j_str = atom_type_j[0].label

            grps.append((grp, grpc,
                         basename + '_Sp-' + str(i + 1) + atom_type_i_str + bdict[bd] + str(j + 1) + atom_type_j_str,
                         'bondExt', (i, j)))

        return grps

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecular group.
        """
        cython.declare(atom=GroupAtom)
        for atom in self.vertices:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return ``True`` if the group contains an atom with the label
        `label` and ``False`` otherwise.
        """
        cython.declare(atom=GroupAtom)
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atom in the group that is labeled with the given `label`.
        Raises :class:`ValueError` if no atom in the group has that label.
        """
        cython.declare(atom=GroupAtom, alist=list)
        alist = [atom for atom in self.vertices if atom.label == label]
        if alist == []:
            raise ValueError('No atom in the functional group \n{1}\n has the label '
                             '"{0}".'.format(label,self.toAdjacencyList()))
        return alist

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        cython.declare(atom=GroupAtom)
        labeled = {}
        for atom in self.vertices:
            if atom.label != '':
                if atom.label in labeled:
                    if isinstance(labeled[atom.label], list):
                        labeled[atom.label].append(atom)
                    else:
                        labeled[atom.label] = [labeled[atom.label]]
                        labeled[atom.label].append(atom)
                else:
                    labeled[atom.label] = atom
        return labeled

    def get_element_count(self):
        """
        Returns the element count for the molecule as a dictionary.
        Wildcards are not counted as any particular element.
        """
        from rmgpy.molecule.atomtype import allElements

        element_count = {}
        for atom in self.atoms:
            same = True
            match = None
            for atomtype in atom.atomType:
                if match is None:
                    # This is the first type in the list, so check all elements
                    for element in allElements:
                        if atomtype.isSpecificCaseOf(atomTypes[element]):
                            match = element
                            break
                else:
                    # We've already matched one atomtype, now confirm that the rest are the same
                    if not atomtype.isSpecificCaseOf(atomTypes[match]):
                        same = False
                        break
            # If match is None, then the group is not a specific case of any element
            if match is not None and same:
                if match in element_count:
                    element_count[match] += 1
                else:
                    element_count[match] = 1

        return element_count

    def fromAdjacencyList(self, adjlist):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from .adjlist import fromAdjacencyList
        self.vertices, multiplicity = fromAdjacencyList(adjlist, group=True)
        if multiplicity is not None:
            self.multiplicity = multiplicity
        self.update()
        return self

    def toAdjacencyList(self, label=''):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from .adjlist import toAdjacencyList
        return toAdjacencyList(self.vertices, multiplicity=self.multiplicity, label='', group=True)

    def updateFingerprint(self):
        """
        Update the molecular fingerprint used to accelerate the subgraph
        isomorphism checks.
        """
        cython.declare(atom=GroupAtom)

        self.elementCount = self.get_element_count()
        self.radicalCount = 0
        for atom in self.vertices:
            if len(atom.radicalElectrons) >= 1:
                self.radicalCount += atom.radicalElectrons[0]

    def isIsomorphic(self, other, initialMap=None, saveOrder=False, strict=True):
        """
        Returns ``True`` if two graphs are isomorphic and ``False``
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        if not strict:
            raise NotImplementedError('There is currently no implementation of the strict argument for Group objects.')
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.isIsomorphic(self, other, initialMap, saveOrder=saveOrder)

    def findIsomorphism(self, other, initialMap=None, saveOrder=False, strict=True):
        """
        Returns ``True`` if `other` is isomorphic and ``False``
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        if not strict:
            raise NotImplementedError('There is currently no implementation of the strict argument for Group objects.')
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.findIsomorphism(self, other, initialMap, saveOrder=saveOrder)

    def isSubgraphIsomorphic(self, other, initialMap=None, generateInitialMap=False, saveOrder=False):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` if self is more specific than other.
        The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short)
        cython.declare(a=GroupAtom, L=list)
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))

        group = other

        if generateInitialMap:
            keys = []
            atms = []
            initialMap = dict()
            for atom in self.atoms:
                if atom.label and atom.label != '':
                    L = [a for a in other.atoms if a.label == atom.label]
                    if L == []:
                        return False
                    elif len(L) == 1:
                        initialMap[atom] = L[0]
                    else:
                        keys.append(atom)
                        atms.append(L)
            if atms:
                for atmlist in itertools.product(*atms):
                    if len(set(atmlist)) != len(atmlist):
                        # skip entries that map multiple graph atoms to the same subgraph atom
                        continue
                    for i, key in enumerate(keys):
                        initialMap[key] = atmlist[i]
                    if (self.isMappingValid(other, initialMap, equivalent=False) and
                            Graph.isSubgraphIsomorphic(self, other, initialMap, saveOrder=saveOrder)):
                        return True
                else:
                    return False
            else:
                if not self.isMappingValid(other, initialMap, equivalent=False):
                    return False

        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return False
        else:
            if group.multiplicity: return False
        # Do the isomorphism comparison
        return Graph.isSubgraphIsomorphic(self, other, initialMap, saveOrder=saveOrder)

    def findSubgraphIsomorphisms(self, other, initialMap=None, saveOrder=False):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` is self is more specific than other.
        Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short)

        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other

        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return []
        else:
            if group.multiplicity:
                return []

        # Do the isomorphism comparison
        return Graph.findSubgraphIsomorphisms(self, other, initialMap, saveOrder=saveOrder)

    def isIdentical(self, other, saveOrder=False):
        """
        Returns ``True`` if `other` is identical and ``False`` otherwise.
        The function `isIsomorphic` respects wildcards, while this function
        does not, make it more useful for checking groups to groups (as
        opposed to molecules to groups)
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # An identical group is always a child of itself and 
        # is the only case where that is true. Therefore
        # if we do both directions of isSubgraphIsmorphic, we need
        # to get True twice for it to be identical
        if not self.isSubgraphIsomorphic(other, None, saveOrder=saveOrder):
            return False
        elif not other.isSubgraphIsomorphic(self, None, saveOrder=saveOrder):
            return False
        else:
            return True

    def isAromaticRing(self):
        """
        This method returns a boolean telling if the group has a 5 or 6 cyclic with
        benzene bonds exclusively
        """

        ring_size = len(self.atoms)
        if ring_size not in [5, 6]:
            return False
        for ringAtom in self.atoms:
            for bondedAtom, bond in ringAtom.edges.items():
                if bondedAtom in self.atoms:
                    if not bond.isBenzene():
                        return False
        return True

    def standardizeAtomType(self):
        """
        This function changes the atomTypes in a group if the atom must
        be a specific atomType based on its bonds and valency.

        Currently only standardizes oxygen, carbon and sulfur atomTypes

        We also only check when there is exactly one atomType,
        one bondType, one radical setting.
        For any group where there are wildcards or multiple attributes,
        we cannot apply this check.

        In the case where the atomType is ambigious based on bonds
        and valency, this function will not change the type.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """
        modified = False

        # If this atom or any of its ligands has wild cards, then don't try to standardize
        if self.hasWildCards: return modified
        for bond12, atom2 in self.bonds.items():
            if atom2.hasWildCards: return modified

        # list of :class:AtomType which are elements with more sub-divided atomtypes beneath them
        specifics = [elementLabel for elementLabel in allElements if elementLabel not in nonSpecifics]
        for atom in self.atoms:
            claimed_atom_type = atom.atomType[0]
            new_atom_type = None
            element = None
            # Ignore elements that do not have more than one atomtype
            if claimed_atom_type.label in nonSpecifics: continue
            for elementLabel in specifics:
                if (claimed_atom_type.label == elementLabel or
                        atomTypes[claimed_atom_type.label] in atomTypes[elementLabel].specific):
                    element = atomTypes[elementLabel]
                    break

            # claimed_atom_type is not in one of the specified elements
            if not element:
                continue
            # Don't standardize atomtypes for nitrogen for now
            # The work on the nitrogen atomtypes is still incomplete
            elif element is atomTypes['N']:
                continue

            group_features = getFeatures(atom, atom.bonds)

            bond_order = atom.getBondOrdersForAtom()
            filled_valency = atom.radicalElectrons[0] + bond_order

            # For an atomtype to be known for certain, the valency must be filled
            # within 1 of the total valency available
            if filled_valency >= PeriodicSystem.valence_electrons[self.symbol][element] - 1:
                for specific_atom_type in element.specific:
                    atomtype_feature_list = specific_atom_type.getFeatures()
                    for mol_feature, atomtype_feature in zip(group_features, atomtype_feature_list):
                        if atomtype_feature == []:
                            continue
                        elif mol_feature not in atomtype_feature:
                            break
                    else:
                        if specific_atom_type is atomTypes['Oa'] or specific_atom_type is atomTypes['Sa']:
                            if atom.lonePairs == 3 or atom.radicalElectrons == 2:
                                new_atom_type = specific_atom_type
                                break
                        else:
                            new_atom_type = specific_atom_type
                            break

            # set the new atom type if the algorithm found one
            if new_atom_type and new_atom_type is not claimed_atom_type:
                atom.atomType[0] = new_atom_type
                modified = True

        return modified

    def createAndConnectAtom(self, atomtypes, connectingAtom, bondOrders):
        """
        This method creates an non-radical, uncharged, :class:GroupAtom with specified list of atomtypes and
        connects it to one atom of the group, 'connectingAtom'. This is useful for making sample atoms.

        Args:
            atomtypes: list of atomtype labels (strs)
            connectingAtom: :class:GroupAtom that is connected to the new benzene atom
            bondOrders: list of bond Orders connecting new_atom and connectingAtom

        Returns: the newly created atom
        """
        atomtypes = [atomTypes[label] for label in atomtypes]  # turn into :class: atomtype instead of labels

        new_atom = GroupAtom(atomType=atomtypes, radicalElectrons=[0], charge=[], label='', lonePairs=None)
        new_bond = GroupBond(connectingAtom, new_atom, order=bondOrders)
        self.addAtom(new_atom)
        self.addBond(new_bond)
        return new_atom

    def addExplicitLigands(self):
        """
        This function O2d/S2d ligand to CO or CS atomtypes if they are not already there.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        atoms_to_add_to = []

        for index, atom in enumerate(self.atoms):
            claimed_atom_type = atom.atomType[0]
            # Do not perform is this atom has wildCards
            if atom.hasWildCards:
                continue
            elif claimed_atom_type is atomTypes['CO'] or claimed_atom_type is atomTypes['CS']:
                for bond12 in atom.bonds.values():
                    if bond12.isDouble():
                        break
                else:
                    atoms_to_add_to.append(index)

        for atomIndex in atoms_to_add_to:
            modified = True
            atomtypes = None
            if self.atoms[atomIndex].atomType[0] is atomTypes['CO']:
                atomtypes = ['O2d']
            elif self.atoms[atomIndex].atomType[0] is atomTypes['CS']:
                atomtypes = ['S2d']
            self.createAndConnectAtom(atomtypes, self.atoms[atomIndex], [2])

        return modified

    def standardizeGroup(self):
        """
        This function modifies groups to make them have a standard AdjList form.

        Currently it makes atomtypes as specific as possible and makes CO/CS atomtypes
        have explicit O2d/S2d ligands. Other functions can be added as necessary

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        # If viable then we apply current conventions:
        check_list = []
        check_list.append(self.standardizeAtomType())
        check_list.append(self.addExplicitLigands())
        return any(check_list)

    def addImplicitAtomsFromAtomType(self):
        """

        Returns: a modified group with implicit atoms added
        Add implicit double/triple bonded atoms O, S or R, for which we will use a C

        Not designed to work with wildcards
        """

        # dictionary of implicit atoms and their bonds
        implicit_atoms = {}
        lone_pairs_required = {}

        copy_group = deepcopy(self)

        for atom1 in copy_group.atoms:
            atomtype_feature_list = atom1.atomType[0].getFeatures()
            lone_pairs_required[atom1] = atomtype_feature_list[8]

            # set to 0 required if empty list
            atomtype_feature_list = [featureList if featureList else [0] for featureList in atomtype_feature_list]
            all_double_required = atomtype_feature_list[1]
            r_double_required = atomtype_feature_list[2]
            o_double_required = atomtype_feature_list[3]
            s_double_required = atomtype_feature_list[4]
            triple_required = atomtype_feature_list[5]
            quadruple_required = atomtype_feature_list[6]

            # count up number of bonds
            single = r_double = o_double = s_double = triple = quadruple = benzene = 0  # note that 0 is immutable
            for atom2, bond12 in atom1.bonds.items():
                # Count numbers of each higher-order bond type
                if bond12.isSingle():
                    single += 1
                elif bond12.isDouble():
                    if atom2.isOxygen():
                        o_double += 1
                    elif atom2.isSulfur():
                        s_double += 1
                    else:
                        # r_double is for double bonds NOT to oxygen or Sulfur
                        r_double += 1
                elif bond12.isTriple():
                    triple += 1
                elif bond12.isQuadruple():
                    quadruple += 1
                elif bond12.isBenzene():
                    benzene += 1

            while o_double < o_double_required[0]:
                o_double += 1
                new_atom = GroupAtom(atomType=[atomTypes['O']], radicalElectrons=[0], charge=[], label='',
                                     lonePairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[2])
                implicit_atoms[new_atom] = new_bond
            while s_double < s_double_required[0]:
                s_double += 1
                new_atom = GroupAtom(atomType=[atomTypes['S']], radicalElectrons=[0], charge=[], label='',
                                     lonePairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[2])
                implicit_atoms[new_atom] = new_bond
            while r_double < r_double_required[0] or r_double + o_double + s_double < all_double_required[0]:
                r_double += 1
                new_atom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='',
                                     lonePairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[2])
                implicit_atoms[new_atom] = new_bond
            while triple < triple_required[0]:
                triple += 1
                new_atom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='',
                                     lonePairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[3])
                implicit_atoms[new_atom] = new_bond
            while quadruple < quadruple_required[0]:
                quadruple += 1
                new_atom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='',
                                     lonePairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[4])
                implicit_atoms[new_atom] = new_bond

        for atom, bond in implicit_atoms.items():
            copy_group.addAtom(atom)
            copy_group.addBond(bond)

        for atom, lonePair in lone_pairs_required.items():
            if lonePair: atom.lonePairs = lonePair

        return copy_group

    def classifyBenzeneCarbons(self, partners=None):
        """
        Args:
            group: :class:Group with atoms to classify
            partners: dictionary of partnered up atoms, which must be a cbf atom

        Returns: tuple with lists of each atom classification
        """
        if not partners:
            partners = {}

        cb_atom_list = []
        cbf_atom_list = []  # All Cbf Atoms
        cbf_atom_list1 = []  # Cbf Atoms that are bonded to exactly one other Cbf (part of 2 rings)
        cbf_atom_list2 = []  # Cbf that are sandwiched between two other Cbf (part of 2 rings)
        connected_cbfs = {}  # dictionary of connections to other cbfAtoms

        # Only want to work with benzene bonds on carbon
        labels_of_carbon_atom_types = [x.label for x in atomTypes['C'].specific] + ['C']
        # Also allow with R!H and some nitrogen groups
        labels_of_carbon_atom_types.extend(['R!H', 'N5b', 'N3b'])

        for atom in self.atoms:
            if atom.atomType[0].label not in labels_of_carbon_atom_types:
                continue
            elif atom.atomType[0].label in ['Cb', 'N5b', 'N3b']:  # Make Cb and N3b into normal cb atoms
                cb_atom_list.append(atom)
            elif atom.atomType[0].label == 'Cbf':
                cbf_atom_list.append(atom)
            else:
                benzene_bonds = 0
                for atom2, bond12 in atom.bonds.items():
                    if bond12.isBenzene(): benzene_bonds += 1
                if benzene_bonds > 2:
                    cbf_atom_list.append(atom)
                elif benzene_bonds > 0:
                    cb_atom_list.append(atom)

        # further sort the cbf atoms
        for cbf_atom in cbf_atom_list:
            fb_bonds = 0
            connected_cbfs[cbf_atom] = []
            for atom2, bond in cbf_atom.bonds.items():
                if bond.order[0] == 1.5 and atom2 in cbf_atom_list:
                    fb_bonds += 1
                    connected_cbfs[cbf_atom].append(atom2)
            if fb_bonds < 2:
                cbf_atom_list1.append(cbf_atom)
            elif fb_bonds == 2:
                cbf_atom_list2.append(cbf_atom)
            elif fb_bonds == 3:
                pass  # leaving here in case we ever want to handle Cbf3 atoms

        # reclassify any atoms with partners as cbf1 atoms
        for cbf_atom in partners:
            if cbf_atom in cb_atom_list:
                cb_atom_list.remove(cbf_atom)
                cbf_atom_list.append(cbf_atom)
                cbf_atom_list1.append(cbf_atom)

        # check that cbfAtoms only have benzene bonds
        for cbf_atom in cbf_atom_list:
            for atom2, bond12 in cbf_atom.bonds.items():
                assert bond12.isBenzene(), "Cbf atom in {0} has a bond with an order other than 1.5".format(self)

        return cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs

    def addImplicitBenzene(self):
        """
        Returns: A modified group with any implicit benzene rings added

        This method currently does not if there are wildcards in atomtypes or bond orders
        The current algorithm also requires that all Cb and Cbf are atomtyped

        There are other cases where the algorithm doesn't work. For example whenever there
        are many dangling Cb or Cbf atoms not in a ring, it is likely fail. In the database test
        (the only use thus far), we will require that any group with more than 3 Cbfs have
        complete rings. This is much stricter than this method can handle, but right now
        this method cannot handle very general cases, so it is better to be conservative. 
        """

        # First define some helper functions
        def checkSet(super_list, sub_list):
            """
            Args:
                super_list: list to check if superset of partList
                sub_list:  list to check if subset of superList

            Returns: Boolean to see if super_list is a superset of sub_list

            """
            super_set = set(super_list)
            sub_set = set(sub_list)
            return super_set.issuperset(sub_set)

        def addCbAtomToRing(ring, cb_atom):
            """
            Every 'Cb' atom belongs in exactly one benzene ring. This function checks
            adds the cb_atom to the ring (in connectivity order) if the cb_atom is connected
            to any the last or first atom in the partial ring.

            Args:
                ring: list of :class:GroupAtoms representing a partial ring to merge
                cb_atom: :class:GroupAtom with atomtype 'Cb'

            Returns: If cb_atom connects to the beginning or end of ring, returns a
            new list of the merged ring, otherwise an empty list

            """

            merged_ring = []
            # ring already complete
            if len(ring) == 6: return merged_ring
            for atom2, bond12 in cb_atom.bonds.items():
                if bond12.isBenzene():
                    if atom2 is ring[-1]:
                        merged_ring = ring + [cb_atom]
                    elif atom2 is ring[0]:
                        merged_ring = [cb_atom] + ring

            return merged_ring

        def mergeOverlappingBenzeneRings(ring1, ring2, od):
            """
            The input arguements of rings are always in the order that the atoms appear
            inside the ring. That is, each atom is connected to the ones adjacent on the
            list.

            Args:
                ring1: list of :class:GroupAtoms representing first partial ring to merge
                ring2: list :class:GroupAtoms representing second partial ring to merge
                od: in for overlap distance

            This function tries to see if the beginning or ends of each list have the
            same atom objects, i.e the two part rings should be merged together.

            Returns: If rings are mergable, returns a new list of the merged ring, otherwise
            an empty list

            """
            new_ring = []
            # ring already complete
            if len(ring1) == 6 or len(ring2) == 6: return new_ring

            # start of ring1 matches end of ring2
            match_list1 = [x1 is x2 for x1, x2 in zip(ring1[-od:], ring2[:od])]
            # end of ring1 matches end of ring2
            match_list2 = [x1 is x2 for x1, x2 in zip(ring1[-od:], ring2[:od - 1:-1])]
            # start of ring1 matches end of ring2
            match_list3 = [x1 is x2 for x1, x2 in zip(ring1[:od], ring2[-od:])]
            # start of ring1 matches start of ring2
            match_list4 = [x1 is x2 for x1, x2 in zip(ring1[:od], ring2[od::-1])]
            if False not in match_list1:
                new_ring = ring1 + ring2[od:]
            elif False not in match_list2:
                new_ring = ring1 + ring2[-od - 1::-1]
            elif False not in match_list3:
                new_ring = ring2[:-od] + ring1
            elif False not in match_list4:
                new_ring = ring2[:od - 1:-1] + ring1

            return new_ring

        #######################################################################################
        # start of main algorithm
        copy_group = deepcopy(self)
        """
        Step 1. Classify all atoms as Cb, Cbf1, Cbf2, Cbf3, ignoring all non-benzene carbons

        Every carbon atom in a benzene ring can be defined as one of the following:
        Cb - benzene carbon in exclusively one ring. Can have one single bond
        Cbf - general classification for any benzene carbon that connects two fused benzene rings
        Cbf1 - Cbf that is bonded to exactly one other Cbf, exclusively in two different benzene rings
        Cbf2 - Cbf that is bonded to exactly two other Cbfs, exclusively in two different benzene ring
        Cbf3 - Cbf that is bonded to exactly three other Cbfs, exclusively in three different benzene rings

        The dictionary connected_cbfs has a cbf atom as key and the other Cbf atoms it is connected to as values

        Currently we only allow 3 Cbf atoms, so Cbf3 is not possible.
        """
        (cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs) = copy_group.classifyBenzeneCarbons()

        """
        #Step 2. Partner up each Cbf1 and Cbf2 atom

        For any fused benzene rings, there will always be exactly two Cbf atoms that join the benzne
        rings. Therefore, we can say that every Cbf1 atom has one 'partner' Cbf atom in which it
        share lies in two benzene rings with. If you try to draw a couple example PAHs, you'll
        find that Cbf2 atoms also have one exclusive 'partner', while Cbf3 atoms are actually
        'partnered' with every atom it is bonded to.

        Because of the current restriction of no more than three Cbf atoms, partnering up every
        Cbf1 is sufficient. If we ever decide to relax this restriction, we will need more code
        to partner up Cbf2 and Cbf3 atoms.
        """
        partners = {}  # dictionary of exclusive partners, has 1:2 and 2:1
        for cbf_atom in cbf_atom_list1:
            if cbf_atom in partners:
                continue
            # if cbf_atom has a connected cbf it must be the partner
            elif connected_cbfs[cbf_atom] and connected_cbfs[cbf_atom][0] not in partners:
                partners[cbf_atom] = connected_cbfs[cbf_atom][0]
                partners[connected_cbfs[cbf_atom][0]] = cbf_atom
            else:
                # search for a potential partner out of atoms benzene bonded atoms
                potential_partner = None
                for atom2, bond12 in cbf_atom.bonds.items():
                    if atom2 in partners:
                        continue
                    # Potential partner must not have any bonds except benzene bonds
                    elif bond12.isBenzene():
                        bonds_are_benzene = [True if bond23.isBenzene() else False for bond23 in atom2.bonds.values()]
                        if all(bonds_are_benzene) and 0 in atom2.radicalElectrons:
                            potential_partner = atom2
                # Make a Cb atom the partner, now marking it as a Cbfatom
                if potential_partner:
                    partners[cbf_atom] = potential_partner
                    partners[potential_partner] = cbf_atom
                # otherwise create a new atom to be the partner
                else:
                    new_atom = copy_group.createAndConnectAtom(['Cbf'], cbf_atom, [1.5])
                    partners[cbf_atom] = new_atom
                    partners[new_atom] = cbf_atom

        # reclassify all atoms since we may have added new ones
        cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs = copy_group.classifyBenzeneCarbons(partners)

        """
        Step 3. Sort all lists by connectivity

        In the coming steps, we will sort Cb/Cbf atom into their benzene rings. If we cannot
        find a ring to sort an atom into, we will create a new ring containing that atom.
        It is important that we always check atoms that are already connected to existing rings
        before completely disconnected atoms. Otherwise, we will erroneously create new rings.
        """
        cb_atom_list = copy_group.sortByConnectivity(cb_atom_list)
        cbf_atom_list1 = copy_group.sortByConnectivity(cbf_atom_list1)
        cbf_atom_list2 = copy_group.sortByConnectivity(cbf_atom_list2)

        """
        Step 4. Initalize the list of rings with any benzene rings that are already explicitly stated

        The variable rings is a list of lists. Each list in rings represents one full benzene rings,
        so it will eventually have six benzene carbons in it. Each ring's list will have the atoms
        sorted by connectivity, such that any atom is bonded to the atoms preceding and following it
        in the list. The first and last atom of the list will also be bonded together.

        """
        rings = [cycle for cycle in copy_group.getAllCyclesOfSize(6) if Group(atoms=cycle).isAromaticRing()]

        """
        Step 5. Add Cbf2 atoms to the correct rings

        In this step, we define 'ring seeds', three atom combination unique to a each benzne ring. We
        then try to merge each ring seed into existing rings. If nosuittable ring is found, we create
        a new ring from the ring seed.

        Every Cbf2 atom is in two rings with unique ring seeds defined by partneredCbf-Cbf2atom-other_cbf
        and partneredCbf-Cbf2Atom-CbAtom. We may need to create the Cbatom in the last ring seed if it
        is not available.
        """
        for cbf_atom in cbf_atom_list2:
            if connected_cbfs[cbf_atom][0] is partners[cbf_atom]:
                other_cbf = connected_cbfs[cbf_atom][1]
            else:
                other_cbf = connected_cbfs[cbf_atom][0]
            # These two ring seeds represent the two unique rings
            new_ring_seeds = [[partners[cbf_atom], cbf_atom, other_cbf],
                              [partners[cbf_atom], cbf_atom]]
            all_ligands = list(cbf_atom.bonds)
            # add a new cb atom to the second new_ring seed
            if len(all_ligands) == 2:
                new_atom = copy_group.createAndConnectAtom(['Cb'], cbf_atom, [1.5])
                new_ring_seeds[1].append(new_atom)
            # join the existing atom to the ringSeed
            elif len(all_ligands) == 3:
                for atom2 in all_ligands:
                    if atom2 not in connected_cbfs[cbf_atom]:
                        new_ring_seeds[1].append(atom2)
                        break
            # Check for duplicates, merge or create new rings
            for index1, ring1 in enumerate(new_ring_seeds):
                merge_ring_dict = {}
                for index, ring2 in enumerate(rings):
                    # check if a duplicate of a fully created ring
                    if checkSet(ring2, ring1): break
                    # Next try to merge the ringseed into rings
                    merge_ring = mergeOverlappingBenzeneRings(ring2, ring1, 2)
                    if merge_ring:
                        merge_ring_dict[index] = merge_ring
                        break
                # otherwise add this ringSeed because it represents a completely new ring
                else:
                    rings.append(ring1)
                # if we merged a ring, we need to remove the old ring from rings and add the merged ring
                rings = [rings[index] if index not in merge_ring_dict else merge_ring_dict[index]
                         for index in range(len(rings))]

        """
        Step 6. Add Cbf1 atoms to the correct rings

        Every Cbf1 atom is in two rings with its partner. In this step, we add this ring seed
        twice to the rings.
        """
        for cbf_atom in cbf_atom_list1:
            new_ring_seed = [partners[cbf_atom], cbf_atom]
            in_ring = 0
            # check to see if duplicate of an existing ring
            for ring in rings:
                if checkSet(ring, new_ring_seed):
                    in_ring += 1
            # move on to next cbf_atom if we found two rings
            if in_ring == 2: continue
            # try to merge into existing rings, if cbf1 is connected
            for index, ring in enumerate(rings):
                merge_ring_dict = {}
                merge_ring = mergeOverlappingBenzeneRings(ring, new_ring_seed, 1)
                if merge_ring:
                    in_ring += 1
                    merge_ring_dict[index] = merge_ring
                    # move on to next cbf_atom if we found two rings
                    if in_ring == 2: break
            # if we merged a ring, we need to remove the old ring from rings and add the merged ring
            rings = [rings[index] if index not in merge_ring_dict else merge_ring_dict[index] for index in
                     range(len(rings))]
            # if we still dont have two ring, we create a completely new ring
            if in_ring < 2:
                for x in range(2 - in_ring):
                    rings.append(copy(new_ring_seed))

        """
        Step 7. Add Cb atoms to the correct rings

        Each Cb atom is part of exactly one benzene ring. In this step we merge or make
        a new ring for each Cb atom.
        """
        for cb_atom in cb_atom_list:
            in_ring = 0
            # check to see if already in a ring
            for ring in rings:
                if checkSet(ring, [cb_atom]):
                    in_ring += 1
            # move on to next ring cb_atom if in a ring
            if in_ring == 1: continue
            # check to see if can be merged to an existing ring
            for index, ring in enumerate(rings):
                merge_ring_dict = {}
                merge_ring = addCbAtomToRing(ring, cb_atom)
                if merge_ring:
                    in_ring += 1
                    merge_ring_dict[index] = merge_ring
                    break
            # if we merged a ring, we need to remove the old ring from rings and add the merged ring
            rings = [rings[index] if index not in merge_ring_dict else merge_ring_dict[index] for index in
                     range(len(rings))]
            # Start completely new ring if not any of above true
            if in_ring == 0:
                rings.append([cb_atom])

        """
        Step 8. Grow each partial ring up to six carbon atoms

        In this step we create new Cb atoms and add them to any rings which do not have 6 atoms.
        """
        merged_ring_dict = {}
        for index, ring in enumerate(rings):
            carbons_to_grow = 6 - len(ring)
            merged_ring_dict[index] = []
            for x in range(carbons_to_grow):
                if x == 0:
                    last_atom = ring[-1]
                else:
                    last_atom = merged_ring_dict[index][-1]
                # add a new atom to the ring and the group
                new_atom = copy_group.createAndConnectAtom(['Cb'], last_atom, [1.5])
                merged_ring_dict[index].append(new_atom)
                # At the end attach to the other endpoint
                if x == carbons_to_grow - 1:
                    new_bond = GroupBond(ring[0], new_atom, order=[1.5])
                    copy_group.addBond(new_bond)

        return copy_group

    def pickWildcards(self):
        """
        Returns: the :class:Group object without wildcards in either atomtype or bonding

        This function will naively pick the first atomtype for each atom, but will try
        to pick bond orders that make sense given the selected atomtypes
        """
        for atom1 in self.atoms:
            atom1.atomType = [atom1.atomType[0]]
            for atom2, bond12 in atom1.bonds.items():
                # skip dynamic bond ordering if there are no wildcards
                if len(bond12.order) < 2:
                    continue
                atom1_features = atom1.atomType[0].getFeatures()
                atom2_features = atom2.atomType[0].getFeatures()
                # skip dynamic bond ordering if there are no features required by the atomtype
                if not any(atom1_features) and not any(atom2_features):
                    bond12.order = [bond12.order[0]]
                    atom2.bonds[atom1].order = bond12.order
                    continue

                atom1_bonds = atom1.countBonds()  # countBonds list must match getFeatures list
                atom2_bonds = atom2.countBonds()
                required_features1 = [atom1_features[x][0] - atom1_bonds[x] if atom1_features[x] else 0 for x in
                                      range(len(atom1_bonds))]
                required_features2 = [atom2_features[x][0] - atom2_bonds[x] if atom2_features[x] else 0 for x in
                                      range(len(atom2_bonds))]

                # subtract 1 from allDouble for each sDouble, oDouble, rDouble so that we don't count them twice
                required_features1[1] = required_features1[1] - required_features1[2] - required_features1[3] - \
                                        required_features1[4]
                required_features2[1] = required_features2[1] - required_features2[2] - required_features2[3] - \
                                        required_features2[4]
                # reverse it because coincidentally the reverse has good priority on what to check first
                required_features1.reverse()
                required_features2.reverse()

                # required features are a now list of [benzene, quadruple, triple, sDouble, oDouble, rDouble, allDouble, single]
                for index, (feature1, feature2) in enumerate(zip(required_features1[:-1], required_features2[:-1])):
                    if feature1 > 0 or feature2 > 0:
                        if index == 0 and 1.5 in bond12.order:  # benzene bonds
                            bond12.order = [1.5]
                            atom2.bonds[atom1].order = bond12.order
                            break
                        elif index == 1 and 4 in bond12.order:  # quadruple bond
                            bond12.order = [4]
                            atom2.bonds[atom1].order = bond12.order
                            break
                        elif index == 2 and 3 in bond12.order:  # triple bond
                            bond12.order = [3]
                            atom2.bonds[atom1].order = bond12.order
                            break
                        elif index > 2 and 2 in bond12.order:  # any case of double bonds
                            if index == 3:  # sDouble bonds
                                if (feature1 > 0 and atom2.isSulfur()) or (feature2 > 0 and atom1.isSulfur()):
                                    bond12.order = [2]
                                    atom2.bonds[atom1].order = bond12.order
                                    break
                            elif index == 4:  # oDoubleBonds
                                if (feature1 > 0 and atom2.isOxygen()) or (feature2 > 0 and atom1.isOxygen()):
                                    bond12.order = [2]
                                    atom2.bonds[atom1].order = bond12.order
                                    break
                            else:  # rDouble or allDouble necessary
                                bond12.order = [2]
                                atom2.bonds[atom1].order = bond12.order
                                break
                else:  # no features required, then pick the first order
                    bond12.order = [bond12.order[0]]
                    atom2.bonds[atom1].order = bond12.order

        # if we have wildcard atomtypes pick one based on ordering of allElements
        for atom in self.atoms:
            for elementLabel in allElements:
                if atomTypes[elementLabel] in atom.atomType[0].specific:
                    atom.atomType = [atomTypes[elementLabel]]
                    break

    def makeSampleMolecule(self):
        """
        Returns: A sample class :Molecule: from the group
        """

        modified_group = self.copy(deep=True)

        # Remove all wildcards
        modified_group.pickWildcards()

        # check that there are less than three Cbf atoms
        cbf_count = 0
        for atom in modified_group.atoms:
            if atom.atomType[0] is atomTypes['Cbf']: cbf_count += 1
        if cbf_count > 3:
            if not modified_group.isBenzeneExplicit():
                raise ImplicitBenzeneError("{0} has more than three Cbf atoms and does not have fully explicit "
                                           "benzene rings.")

        # Add implicit atoms
        modified_group = modified_group.addImplicitAtomsFromAtomType()

        # Add implicit benzene rings
        if not modified_group.isBenzeneExplicit():
            modified_group = modified_group.addImplicitBenzene()
        # Make dictionary of :GroupAtoms: to :Atoms: and vice versa
        group_to_mol = {}
        mol_to_group = {}
        for atom in modified_group.atoms:
            mol_atom = atom.makeSampleAtom()
            group_to_mol[atom] = mol_atom
            mol_to_group[mol_atom] = atom

        # create the molecule
        new_molecule = mol.Molecule(atoms=list(group_to_mol.values()))

        # Add explicit bonds to :Atoms:
        for atom1 in modified_group.atoms:
            for atom2, bond12 in atom1.bonds.items():
                bond12.makeBond(new_molecule, group_to_mol[atom1], group_to_mol[atom2])

        # Saturate up to expected valency
        for mol_atom in new_molecule.atoms:
            if mol_atom.charge:
                stated_charge = mol_atom.charge
            # otherwise assume no charge (or implicit atoms we assume hvae no charge)
            else:
                stated_charge = 0
            mol_atom.updateCharge()
            if mol_atom.charge - stated_charge:
                hydrogen_needed = mol_atom.charge - stated_charge
                if mol_atom in mol_to_group and mol_to_group[mol_atom].atomType[0].single:
                    max_single = max(mol_to_group[mol_atom].atomType[0].single)
                    single_present = sum([1 for atom in mol_atom.bonds if mol_atom.bonds[atom].isSingle()])
                    max_hydrogen = max_single - single_present
                    if hydrogen_needed > max_hydrogen: hydrogen_needed = max_hydrogen
                for x in range(hydrogen_needed):
                    new_h = mol.Atom('H', radicalElectrons=0, lonePairs=0, charge=0)
                    new_bond = mol.Bond(mol_atom, new_h, 1)
                    new_molecule.addAtom(new_h)
                    new_molecule.addBond(new_bond)
                mol_atom.updateCharge()

        new_molecule.update()

        # Check that the charge of atoms is expected
        for atom in new_molecule.atoms:
            if abs(atom.charge) > 0:
                if atom in mol_to_group:
                    group_atom = mol_to_group[atom]
                else:
                    raise UnexpectedChargeError(graph=new_molecule)
                # check hardcoded atomtypes
                positive_charged = ['Csc', 'Cdc',
                                    'N3sc', 'N5sc', 'N5dc', 'N5ddc', 'N5tc', 'N5b',
                                    'O2sc', 'O4sc', 'O4dc', 'O4tc',
                                    'S2sc', 'S4sc', 'S4dc', 'S4tdc', 'S6sc', 'S6dc', 'S6tdc']
                negative_charged = ['C2sc', 'C2dc', 'C2tc',
                                    'N0sc', 'N1sc', 'N1dc', 'N5dddc',
                                    'O0sc',
                                    'S0sc', 'S2sc', 'S2dc', 'S2tc', 'S4dc', 'S4tdc', 'S6sc', 'S6dc', 'S6tdc']
                if group_atom.atomType[0] in [atomTypes[x] for x in positive_charged] and atom.charge > 0:
                    pass
                elif group_atom.atomType[0] in [atomTypes[x] for x in negative_charged] and atom.charge < 0:
                    pass
                # declared charge in original group is not same as new charge
                elif atom.charge in group_atom.charge:
                    pass
                else:
                    raise UnexpectedChargeError(graph=new_molecule)

        return new_molecule

    def isBenzeneExplicit(self):
        """

        Returns: 'True' if all Cb, Cbf atoms are in completely explicitly stated benzene rings.

        Otherwise return 'False'

        """

        # classify atoms
        cb_atom_list = []

        # only want to work with carbon atoms
        labels_of_carbon_atom_types = [x.label for x in atomTypes['C'].specific] + ['C', 'N3b', 'N5b']

        for atom in self.atoms:
            if atom.atomType[0].label not in labels_of_carbon_atom_types:
                continue
            elif atom.atomType[0].label in ['Cb', 'Cbf', 'N3b', 'N5b']:  # Make Cb and N3b into normal cb atoms
                cb_atom_list.append(atom)
            else:
                benzene_bonds = 0
                for atom2, bond12 in atom.bonds.items():
                    if bond12.isBenzene(): benzene_bonds += 1
                if benzene_bonds > 0: cb_atom_list.append(atom)

        # get all explicit benzene rings
        rings = [cycle for cycle in self.getAllCyclesOfSize(6) if Group(atoms=cycle).isAromaticRing()]

        # test that all benzene atoms are in benzene rings
        for atom in cb_atom_list:
            in_ring = False
            for ring in rings:
                if atom in ring:
                    in_ring = True
            if not in_ring:
                return False
        else:
            return True

    def mergeGroups(self, other, keepIdenticalLabels=False):
        """
        This function takes `other` :class:Group object and returns a merged :class:Group object based
        on overlapping labeled atoms between self and other

        Currently assumes `other` can be merged at the closest labelled atom
        if keepIdenticalLabels=True mergeGroups will not try to merge atoms with the same labels
        """
        labeled1 = self.getLabeledAtoms()
        labeled2 = other.getLabeledAtoms()
        overlapping_labels = [x for x in labeled1 if x in labeled2]

        # dictionary of key = original atoms, value = copy atoms for deep copies of the two groups
        self_dict = self.copyAndMap()
        other_dict = other.copyAndMap()

        # Sort atoms to go into the merged_group
        merged_group_atoms = list(other_dict.values())  # all end atoms with end up in the merged_group
        # only non-overlapping atoms from the backbone will be in the merged_group
        for originalAtom, newAtom in self_dict.items():
            if originalAtom.label not in overlapping_labels:
                merged_group_atoms.append(newAtom)

        merged_group = Group(atoms=merged_group_atoms)

        """
        The following loop will move bonds that are exclusively in the new backbone so that
        they connect to the backbone. For example, assume backbone has bond atomA-atomB,
        where atomB is labelled as *2 and there is an atomC in the end analgously labelled
        *2. We need to remove the bond between atomA and atomB. Then we need to add a bond
        between atomA and atomC.
        """
        if not keepIdenticalLabels:
            bonds_to_remove = []
            for label in overlapping_labels:
                old_atom_b = self.getLabeledAtom(label)[0]
                for old_atom_a, old_bond_ab in old_atom_b.bonds.items():
                    if old_atom_a.label not in overlapping_labels:  # this is bond we need to transfer over
                        # find and record bondAB from new backbone for later removal
                        new_atom_a = self_dict[old_atom_a]
                        new_atom_b = self_dict[old_atom_b]
                        new_atom_c = merged_group.getLabeledAtom(old_atom_b.label)[0]
                        for atom, new_bond_ab in new_atom_a.bonds.items():
                            if atom is new_atom_b:
                                bonds_to_remove.append(new_bond_ab)
                                break
                        # add bond between atomA and AtomC
                        new_bond_ac = GroupBond(new_atom_a, new_atom_c, order=old_bond_ab.order)
                        merged_group.addBond(new_bond_ac)
            # remove bonds from merged_group
            for bond in bonds_to_remove:
                merged_group.removeBond(bond)

        return merged_group

    def resetRingMembership(self):
        """
        Resets ring membership information in the GroupAtom.props attribute.
        """
        cython.declare(atom=GroupAtom)

        for atom in self.atoms:
            if 'inRing' in atom.props:
                del atom.props['inRing']
