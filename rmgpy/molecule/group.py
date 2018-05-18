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
This module provides classes and methods for working with molecular substructure
groups. These enable molecules to be searched for common motifs (e.g.
reaction sites).
"""

import cython

from .graph import Vertex, Edge, Graph
from .atomtype import atomTypes, allElements, nonSpecifics, getFeatures
from .element import PeriodicSystem
import rmgpy.molecule.molecule as mol
import rmgpy.molecule.element as elements
from copy import deepcopy, copy
from rmgpy.exceptions import ActionError, ImplicitBenzeneError, UnexpectedChargeError

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
    def bonds(self): return self.edges

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
                raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: Invalid order "{0}".'.format(order))
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __formBond(self, order):
        """
        Update the atom group as a result of applying a FORM_BOND action,
        where `order` specifies the order of the forming bond, and should be
        1 (since we only allow forming of single bonds).
        """
        if order != 1:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.formBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __breakBond(self, order):
        """
        Update the atom group as a result of applying a BREAK_BOND action,
        where `order` specifies the order of the breaking bond, and should be
        1 (since we only allow breaking of single bonds).
        """
        if order != 1:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Invalid order "{0}".'.format(order))
        atomType = []
        for atom in self.atomType:
            atomType.extend(atom.breakBond)
        if len(atomType) == 0:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: Unknown atom type produced from set "{0}".'.format(self.atomType))
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
            raise ActionError('Unable to update GroupAtom due to GAIN_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        if not self.radicalElectrons:
            radicalElectrons = [1,2,3,4]
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
            raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))

        if not self.radicalElectrons:
            radicalElectrons = [0,1,2,3]
        else:
            for electron in self.radicalElectrons:
                electron = electron - radical
                if electron < 0:
                    raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))
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
            raise ActionError('Unable to update GroupAtom due to GAIN_PAIR action: Unknown atom type produced from set "{0}".'.format(self.atomType))

        #Add a lone pair to a group atom with none
        if not self.lonePairs:
            self.lonePairs = [1,2,3,4] #set to a wildcard of any number greater than 0
        #Add a lone pair to a group atom that already has at least one lone pair
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
            raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: Unknown atom type produced from set "{0}".'.format(self.atomType))

        if not self.lonePairs:
            self.lonePairs = [0,1,2,3] #set to a wildcard of any number fewer than 4
        else:
            for x in self.lonePairs:
                if x - pair < 0:
                    raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: Invalid lone electron pairs set "{0}".'.format(self.lonePairs))
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

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Atom` or an :class:`GroupAtom`
        object. When comparing two :class:`GroupAtom` objects, this function
        respects wildcards, e.g. ``R!H`` is equivalent to ``C``.
        
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the equivalent method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)
        group=other
        
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
                for radical2  in group.radicalElectrons:
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
        #Each charge in self must have an equivalent in other (and vice versa)
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
        group=other
        
        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short, 
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType: # all these must match
            for atomType2 in group.atomType: # can match any of these
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
        #Each charge in self must have an equivalent in other
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

    def isOxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if
        not.
        """
        allOxygens = [atomTypes['O']] + atomTypes['O'].specific
        checkList=[x in allOxygens for x in self.atomType]

        return all(checkList)

    def isSulfur(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if
        not.
        """
        allSulfur = [atomTypes['S']] + atomTypes['S'].specific
        checkList=[x in allSulfur for x in self.atomType]

        return all(checkList)

    def hasWildcards(self):
        """
        Return ``True`` if the atom has wildcards in any of the attributes:
        atomtype, electronpairs, lone pairs, charge, and bond order. Returns
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

    def countBonds(self, wildcards = False):
        """
        Returns: list of the number of bonds currently on the :class:GroupAtom

        If the argument wildcards is turned off then any bonds with multiple
        options for bond orders will not be counted
        """
        #count up number of bonds
        single = 0; rDouble = 0; oDouble = 0; sDouble = 0; triple = 0; benzene = 0
        for atom2, bond12 in self.bonds.iteritems():
            if not wildcards and len(bond12.order) > 1:
                continue
            # Count numbers of each higher-order bond type
            if bond12.isSingle(wildcards = True):
                single += 1
            if bond12.isDouble(wildcards = True):
                if atom2.isOxygen():
                    oDouble += 1
                elif atom2.isSulfur():
                    sDouble += 1
                else:
                    # rDouble is for double bonds NOT to oxygen or Sulfur
                    rDouble += 1
            if bond12.isTriple(wildcards = True): triple += 1
            if bond12.isBenzene(wildcards = True): benzene += 1

        allDouble = rDouble + oDouble + sDouble

        return [single, allDouble, rDouble, oDouble, sDouble, triple, benzene]

    def makeSampleAtom(self):
        """

        Returns: a class :Atom: object analagous to the GroupAtom

        This makes a sample, so it takes the first element when there are multiple options inside of
        self.atomtype, self.radicalElectrons, self.lonePairs, and self.charge

        """

        #Use the first atomtype to determine element, even if there is more than one atomtype
        atomtype = self.atomType[0]
        element = None

        defaultLonePairs={'H': 0,
                          'D': 0,
                          'T': 0,
                          'He':1,
                          'C': 0,
                          'O': 2,
                          'N': 1,
                          'Si':0,
                          'S': 2,
                          'Ne':4,
                          'Cl':3,
                          'Ar':4,
        }

        for elementLabel in allElements:
            if atomtype is atomTypes[elementLabel] or atomtype in atomTypes[elementLabel].specific:
                element = elementLabel
                break
        else:
            #For types that correspond to more than one type of element, pick the first that appears in specific
            for subtype in atomtype.specific:
                if subtype.label in allElements:
                    element = subtype.label
                    break

        #dummy defaultAtom to get default values
        defaultAtom = mol.Atom()

        #Three possible values for charge and lonePairs
        if self.charge:
            newCharge = self.charge[0]
        elif atomtype.charge:
            newCharge = atomtype.charge[0]
        else:
            newCharge = defaultAtom.charge

        if self.lonePairs:
            newLonePairs = self.lonePairs[0]
        elif atomtype.lonePairs:
            newLonePairs = atomtype.lonePairs[0]
        else:
            newLonePairs = defaultAtom.lonePairs

        newAtom = mol.Atom(element = element,
                           radicalElectrons = self.radicalElectrons[0] if self.radicalElectrons else defaultAtom.radicalElectrons,
                           charge = newCharge,
                           lonePairs = newLonePairs,
                           label = self.label if self.label else defaultAtom.label)

        #For some reason the default when no lone pairs is set to -100,
        #Based on git history, it is probably because RDKit requires a number instead of None
        if newAtom.lonePairs == -100:
            newAtom.lonePairs = defaultLonePairs[newAtom.symbol]

        return newAtom

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
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. a bond will match the
    group if it matches *any* item in the list.
    """

    def __init__(self, atom1, atom2, order=None):
        Edge.__init__(self, atom1, atom2)
        if order is not None and all([isinstance(oneOrder,str) for oneOrder in order]):
            self.setOrderStr(order)
        elif order is not None and any([isinstance(oneOrder,str) for oneOrder in order]):
            raise ActionError('order list given {} does not consist of only strings or only numbers'.format(order))
        else:
            self.order = order or []

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
            elif value == 1.5:
                values.append('B')
            elif value == 0:
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
            elif value == 'B':
                values.append(1.5)
            elif value == 'H':
                values.append(0)
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
            
    def isSingle(self, wildcards = False):
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
                if abs(order-1) <= 1e-9:
                    return True
            else: return False
        else:
            return abs(self.order[0]-1) <= 1e-9 and len(self.order) == 1

    def isDouble(self, wildcards = False):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are double.
        """
        if wildcards:
            for order in self.order:
                if abs(order-2) <= 1e-9:
                    return True
            else: return False
        else:
            return abs(self.order[0]-2) <= 1e-9 and len(self.order) == 1

    def isTriple(self, wildcards = False):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are triple.
        """
        if wildcards:
            for order in self.order:
                if abs(order-3) <= 1e-9:
                    return True
            else: return False
        else:
            return abs(self.order[0]-3) <= 1e-9 and len(self.order) == 1

    def isBenzene(self, wildcards = False):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are benzene
        """
        if wildcards:
            for order in self.order:
                if abs(order-1.5) <= 1e-9:
                    return True
            else: return False
        else:
            return abs(self.order[0]-1.5) <= 1e-9 and len(self.order) == 1

    def isHydrogenBond(self, wildcards = False):
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
            else: return False
        else:
            return abs(self.order[0]) <= 1e-9 and len(self.order) == 1
        
    def __changeBond(self, order):
        """
        Update the bond group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order. `order` is normally 1 or -1, but can be any value
        """
        newOrder = [value + order for value in self.order]
        if any([value < 0 or value > 3 for value in newOrder]):
            raise ActionError('Unable to update Bond due to CHANGE_BOND action: Invalid resulting order "{0}".'.format(newOrder))
        # Change any modified benzene orders to the appropriate stable order
        newOrder = set(newOrder)
        if 0.5 in newOrder:
            newOrder.remove(0.5)
            newOrder.add(1)
        if 2.5 in newOrder:
            newOrder.remove(2.5)
            newOrder.add(2)
        # Set the new bond orders
        self.order = list(newOrder)

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
        for order1 in self.order: # all these must match
            for order2 in gb.order: # can match any of these
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
        newBond = mol.Bond(atom1, atom2, order = self.order[0])
        molecule.addBond(newBond)

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
            atomType = '{0!s} '.format(atom.label if atom.label != '' else '')
            atomType += ','.join([at.label for at in atom.atomType])
            atomType = '"' + atomType + '"'
            graph.add_node(pydot.Node(name=str(index + 1), label=atomType, fontname="Helvetica", fontsize="16"))
        for atom1 in self.atoms:
            for atom2, bond in atom1.bonds.iteritems():
                index1 = self.atoms.index(atom1)
                index2 = self.atoms.index(atom2)
                if index1 < index2:
                    bondType = ','.join([order for order in bond.getOrderStr()])
                    bondType = '"' + bondType + '"'
                    graph.add_edge(pydot.Edge(src=str(index1 + 1), dst=str(index2 + 1), label=bondType, fontname="Helvetica", fontsize="16"))

        img = graph.create(prog='neato', format=format)
        return img

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

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

        sortedAtomList=[]
        sortedAtomList.append(atomList.pop(0))
        while atomList:
            for atom1 in sortedAtomList:
                added = False
                for atom2, bond12 in atom1.bonds.iteritems():
                    if bond12.isBenzene() and atom2 in atomList:
                        sortedAtomList.append(atom2)
                        atomList.remove(atom2)
                        added = True
                if added: break
            else:
                sortedAtomList.append(atomList.pop(0))

        return sortedAtomList

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
                        for _ , bond in atom.bonds.iteritems():
                            bond_order += bond.order[0]
                        lone_pairs = atom.lonePairs[0]
                        radical_electrons = atom.radicalElectrons[0]
                        atom.charge[0] = valence_electron - bond_order - 2 * lone_pairs - radical_electrons
                    else:
                        # if the group is not specified to specific element, charge will not be updated
                        pass

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
        cython.declare(atom=GroupAtom)
        for atom in self.vertices:
            if atom.label == label: return atom
        raise ValueError('No atom in the functional group has the label "{0}".'.format(label))

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
                    if isinstance(labeled[atom.label],list):
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

    def isIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if two graphs are isomorphic and ``False``
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.isIsomorphic(self, other, initialMap)

    def findIsomorphism(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is isomorphic and ``False``
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.findIsomorphism(self, other, initialMap)

    def isSubgraphIsomorphic(self, other, initialMap=None):
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
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
        
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
        return Graph.isSubgraphIsomorphic(self, other, initialMap)

    def findSubgraphIsomorphisms(self, other, initialMap=None):
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
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other
        
        if self.multiplicity:
            for mult1 in self.multiplicity:
                if group.multiplicity:
                    for mult2 in group.multiplicity:
                        if mult1 == mult2: break
                    else:
                        return []
        else:
            if group.multiplicity: return []
                
        # Do the isomorphism comparison
        return Graph.findSubgraphIsomorphisms(self, other, initialMap)
    
    def isIdentical(self, other):
        """
        Returns ``True`` if `other` is identical and ``False`` otherwise.
        The function `isIsomorphic` respects wildcards, while this function
        does not, make it more useful for checking groups to groups (as
        opposed to molecules to groups)
        """
        # It only makes sense to compare a Group to a Group for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # An identical group is always a child of itself and 
        # is the only case where that is true. Therefore
        # if we do both directions of isSubgraphIsmorphic, we need
        # to get True twice for it to be identical
        if not self.isSubgraphIsomorphic(other):
            return False
        elif not other.isSubgraphIsomorphic(self):
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
            for bondedAtom, bond in ringAtom.edges.iteritems():
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

        #If this atom or any of its ligands has wild cards, then don't try to standardize
        if self.hasWildCards: return modified
        for bond12, atom2 in self.bonds.iteritems():
            if atom2.hasWildCards: return modified

        #list of :class:AtomType which are elements with more sub-divided atomtypes beneath them
        specifics= [elementLabel for elementLabel in allElements if elementLabel not in nonSpecifics]
        for atom in self.atoms:
            claimedAtomType = atom.atomType[0]
            newAtomType = None
            element = None
            #Ignore elements that do not have more than one atomtype
            if claimedAtomType.label in nonSpecifics: continue
            for elementLabel in specifics:
                if claimedAtomType.label == elementLabel or atomTypes[claimedAtomType.label] in atomTypes[elementLabel].specific:
                    element = atomTypes[elementLabel]
                    break

            #claimedAtomType is not in one of the specified elements
            if not element: continue
            #Don't standardize atomtypes for nitrogen for now
            #The work on the nitrogen atomtypes is still incomplete
            elif element is atomTypes['N']: continue

            groupFeatures = getFeatures(atom, atom.bonds)

            bondOrder = atom.getBondOrdersForAtom()
            filledValency =  atom.radicalElectrons[0] + bondOrder

            #For an atomtype to be known for certain, the valency must be filled
            #within 1 of the total valency available
            if filledValency >= PeriodicSystem.valence_electrons[self.symbol][element] - 1:
                for specificAtomType in element.specific:
                    atomtypeFeatureList = specificAtomType.getFeatures()
                    for molFeature, atomtypeFeature in zip(groupFeatures, atomtypeFeatureList):
                        if atomtypeFeature == []:
                            continue
                        elif molFeature not in atomtypeFeature:
                            break
                    else:
                        if specificAtomType is atomTypes['Oa'] or specificAtomType is atomTypes['Sa']:
                            if atom.lonePairs == 3 or atom.radicalElectrons == 2:
                                newAtomType = specificAtomType
                                break
                        else:
                            newAtomType = specificAtomType
                            break

            #set the new atom type if the algorithm found one
            if newAtomType and newAtomType is not claimedAtomType:
                atom.atomType[0] = newAtomType
                modified = True

        return modified

    def createAndConnectAtom(self, atomtypes, connectingAtom, bondOrders):
        """
        This method creates an non-radical, uncharged, :class:GroupAtom with specified list of atomtypes and
        connects it to one atom of the group, 'connectingAtom'. This is useful for making sample atoms.

        Args:
            atomtypes: list of atomtype labels (strs)
            connectingAtom: :class:GroupAtom that is connected to the new benzene atom
            bondOrders: list of bond Orders connecting newAtom and connectingAtom

        Returns: the newly created atom
        """
        atomtypes = [atomTypes[label] for label in atomtypes] #turn into :class: atomtype instead of labels

        newAtom = GroupAtom(atomType= atomtypes, radicalElectrons=[0], charge=[], label='', lonePairs=None)
        newBond = GroupBond(connectingAtom, newAtom, order=bondOrders)
        self.addAtom(newAtom)
        self.addBond(newBond)
        return newAtom

    def addExplicitLigands(self):
        """
        This function O2d/S2d ligand to CO or CS atomtypes if they are not already there.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        atomsToAddTo=[]

        for index, atom in enumerate(self.atoms):
            claimedAtomType = atom.atomType[0]
            #Do not perform is this atom has wildCards
            if atom.hasWildCards: continue
            elif claimedAtomType is atomTypes['CO'] or claimedAtomType is atomTypes['CS']:
                for bond12 in atom.bonds.itervalues():
                    if bond12.isDouble():
                        break
                else: atomsToAddTo.append(index)

        for atomIndex in atomsToAddTo:
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
        checkList=[]
        checkList.append(self.standardizeAtomType())
        checkList.append(self.addExplicitLigands())
        return any(checkList)


    def addImplicitAtomsFromAtomType(self):
        """

        Returns: a modified group with implicit atoms added
        Add implicit double/triple bonded atoms O, S or R, for which we will use a C

        Not designed to work with wildcards
        """

        #dictionary of implicit atoms and their bonds
        implicitAtoms = {}
        lonePairsRequired = {}

        copyGroup = deepcopy(self)

        for atom1 in copyGroup.atoms:
            atomtypeFeatureList = atom1.atomType[0].getFeatures()
            lonePairsRequired[atom1]=atomtypeFeatureList[7]

            #set to 0 required if empty list
            atomtypeFeatureList = [featureList if featureList else [0] for featureList in atomtypeFeatureList]
            allDoubleRequired = atomtypeFeatureList[1]
            rDoubleRequired = atomtypeFeatureList[2]
            oDoubleRequired = atomtypeFeatureList[3]
            sDoubleRequired = atomtypeFeatureList[4]
            tripleRequired = atomtypeFeatureList[5]

            #count up number of bonds
            single = 0; rDouble = 0; oDouble = 0; sDouble = 0; triple = 0; benzene = 0
            for atom2, bond12 in atom1.bonds.iteritems():
                # Count numbers of each higher-order bond type
                if bond12.isSingle():
                    single += 1
                elif bond12.isDouble():
                    if atom2.isOxygen():
                        oDouble += 1
                    elif atom2.isSulfur():
                        sDouble += 1
                    else:
                        # rDouble is for double bonds NOT to oxygen or Sulfur
                        rDouble += 1
                elif bond12.isTriple(): triple += 1
                elif bond12.isBenzene(): benzene += 1


            while oDouble < oDoubleRequired[0]:
                oDouble +=1
                newAtom = GroupAtom(atomType=[atomTypes['O']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=[2])
                implicitAtoms[newAtom] = newBond
            while sDouble < sDoubleRequired[0]:
                sDouble +=1
                newAtom = GroupAtom(atomType=[atomTypes['S']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=[2])
                implicitAtoms[newAtom] = newBond
            while rDouble < rDoubleRequired[0] or rDouble + oDouble + sDouble < allDoubleRequired[0]:
                rDouble +=1
                newAtom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=[2])
                implicitAtoms[newAtom] = newBond
            while triple < tripleRequired[0]:
                triple +=1
                newAtom = GroupAtom(atomType=[atomTypes['C']], radicalElectrons=[0], charge=[], label='', lonePairs=None)
                newBond = GroupBond(atom1, newAtom, order=[3])
                implicitAtoms[newAtom] = newBond

        for atom, bond in implicitAtoms.iteritems():
            copyGroup.addAtom(atom)
            copyGroup.addBond(bond)

        for atom, lonePair in lonePairsRequired.iteritems():
            if lonePair: atom.lonePairs = lonePair

        return copyGroup

    def classifyBenzeneCarbons(self, partners = None):
        """
        Args:
            group: :class:Group with atoms to classify
            partners: dictionary of partnered up atoms, which must be a cbf atom

        Returns: tuple with lists of each atom classification
        """
        if not partners: partners = {}

        cbAtomList = []
        cbfAtomList = [] #All Cbf Atoms
        cbfAtomList1 = [] #Cbf Atoms that are bonded to exactly one other Cbf (part of 2 rings)
        cbfAtomList2 = [] #Cbf that are sandwiched between two other Cbf (part of 2 rings)
        connectedCbfs={} #dictionary of connections to other cbfAtoms

        #Only want to work with benzene bonds on carbon
        labelsOfCarbonAtomTypes = [x.label for x in atomTypes['C'].specific] + ['C']
        #Also allow with R!H and some nitrogen groups
        labelsOfCarbonAtomTypes.extend(['R!H', 'N5b', 'N3b'])

        for atom in self.atoms:
            if not atom.atomType[0].label in labelsOfCarbonAtomTypes: continue
            elif atom.atomType[0].label in ['Cb', 'N5b', 'N3b']: #Make Cb and N3b into normal cb atoms
                cbAtomList.append(atom)
            elif atom.atomType[0].label == 'Cbf':
                cbfAtomList.append(atom)
            else:
                benzeneBonds = 0
                for atom2, bond12 in atom.bonds.iteritems():
                    if bond12.isBenzene(): benzeneBonds+=1
                if benzeneBonds > 2: cbfAtomList.append(atom)
                elif benzeneBonds >0: cbAtomList.append(atom)

        #further sort the cbf atoms
        for cbfAtom in cbfAtomList:
            fbBonds = 0
            connectedCbfs[cbfAtom] = []
            for atom2, bond in cbfAtom.bonds.iteritems():
                if bond.order[0] == 1.5 and atom2 in cbfAtomList:
                    fbBonds +=1
                    connectedCbfs[cbfAtom].append(atom2)
            if fbBonds < 2: cbfAtomList1.append(cbfAtom)
            elif fbBonds == 2: cbfAtomList2.append(cbfAtom)
            elif fbBonds == 3: pass #leaving here in case we ever want to handle Cbf3 atoms

        #reclassify any atoms with partners as cbf1 atoms
        for cbfAtom in partners:
            if cbfAtom in cbAtomList:
                cbAtomList.remove(cbfAtom)
                cbfAtomList.append(cbfAtom)
                cbfAtomList1.append(cbfAtom)

        #check that cbfAtoms only have benzene bonds
        for cbfAtom in cbfAtomList:
            for atom2, bond12 in cbfAtom.bonds.iteritems():
                assert bond12.isBenzene(), "Cbf atom in {0} has a bond with an order other than 1.5".format(self)

        return (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs)

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
        #First define some helper functions
        def checkSet(superList, subList):
            """
            Args:
                superList: list to check if superset of partList
                subList:  list to check if subset of superList

            Returns: Boolean to see if superList is a superset of subList

            """
            superSet = set(superList)
            subSet = set(subList)
            return superSet.issuperset(subSet)

        def addCbAtomToRing(ring, cbAtom):
            """
            Every 'Cb' atom belongs in exactly one benzene ring. This function checks
            adds the cbAtom to the ring (in connectivity order) if the cbAtom is connected
            to any the last or first atom in the partial ring.

            Args:
                ring: list of :class:GroupAtoms representing a partial ring to merge
                cbAtom: :class:GroupAtom with atomtype 'Cb'

            Returns: If cbAtom connects to the beginning or end of ring, returns a
            new list of the merged ring, otherwise an empty list

            """

            mergedRing = []
            #ring already complete
            if len(ring) == 6 : return mergedRing
            for atom2, bond12 in cbAtom.bonds.iteritems():
                if bond12.isBenzene():
                    if atom2 is ring[-1]: mergedRing = ring+[cbAtom]
                    elif atom2 is ring[0]: mergedRing = [cbAtom] +ring

            return mergedRing

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
            newRing = []
            #ring already complete
            if len(ring1) ==6 or len(ring2) == 6: return newRing

            #start of ring1 matches end of ring2
            matchList1 = [x1 is x2 for x1,x2 in zip(ring1[-od:],ring2[:od])]
            #end of ring1 matches end of ring2
            matchList2 = [x1 is x2 for x1,x2 in zip(ring1[-od:],ring2[:od-1:-1])]
            #start of ring1 matches end of ring2
            matchList3 = [x1 is x2 for x1,x2 in zip(ring1[:od],ring2[-od:])]
            #start of ring1 matches start of ring2
            matchList4 = [x1 is x2 for x1,x2 in zip(ring1[:od],ring2[od::-1])]
            if not False in matchList1:
                newRing = ring1 +ring2[od:]
            elif not False in matchList2:
                newRing = ring1 + ring2[-od-1::-1]
            elif not False in matchList3:
                newRing = ring2[:-od] + ring1
            elif not False in matchList4:
                newRing = ring2[:od-1:-1] + ring1

            return newRing
        #######################################################################################
        #start of main algorithm
        copyGroup = deepcopy(self)
        """
        Step 1. Classify all atoms as Cb, Cbf1, Cbf2, Cbf3, ignoring all non-benzene carbons

        Every carbon atom in a benzene ring can be defined as one of the following:
        Cb - benzene carbon in exclusively one ring. Can have one single bond
        Cbf - general classification for any benzene carbon that connects two fused benzene rings
        Cbf1 - Cbf that is bonded to exactly one other Cbf, exclusively in two different benzene rings
        Cbf2 - Cbf that is bonded to exactly two other Cbfs, exclusively in two different benzene ring
        Cbf3 - Cbf that is bonded to exactly three other Cbfs, exclusively in three different benzene rings

        The dictionary connectedCbfs has a cbf atom as key and the other Cbf atoms it is connected to as values

        Currently we only allow 3 Cbf atoms, so Cbf3 is not possible.
        """
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = copyGroup.classifyBenzeneCarbons()

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
        partners = {} #dictionary of exclusive partners, has 1:2 and 2:1
        for cbfAtom in cbfAtomList1:
            if cbfAtom in partners: continue
            #if cbfAtom has a connected cbf it must be the partner
            elif connectedCbfs[cbfAtom] and connectedCbfs[cbfAtom][0] not in partners:
                partners[cbfAtom] = connectedCbfs[cbfAtom][0]
                partners[connectedCbfs[cbfAtom][0]] = cbfAtom
            else:
                #search for a potential partner out of atoms benzene bonded atoms
                potentialPartner = None
                for atom2, bond12 in cbfAtom.bonds.iteritems():
                    if atom2 in partners: continue
                    #Potential partner must not have any bonds except benzene bonds
                    elif bond12.isBenzene():
                        bondsAreBenzene = [True if bond23.isBenzene() else False for bond23 in atom2.bonds.values()]
                        if all(bondsAreBenzene) and 0 in atom2.radicalElectrons:
                            potentialPartner= atom2
                #Make a Cb atom the partner, now marking it as a Cbfatom
                if potentialPartner:
                    partners[cbfAtom] = potentialPartner
                    partners[potentialPartner]= cbfAtom
                #otherwise create a new atom to be the partner
                else:
                    newAtom = copyGroup.createAndConnectAtom(['Cbf'], cbfAtom, [1.5])
                    partners[cbfAtom] = newAtom
                    partners[newAtom] = cbfAtom

        #reclassify all atoms since we may have added new ones
        (cbAtomList, cbfAtomList, cbfAtomList1, cbfAtomList2, connectedCbfs) = copyGroup.classifyBenzeneCarbons(partners)

        """
        Step 3. Sort all lists by connectivity

        In the coming steps, we will sort Cb/Cbf atom into their benzene rings. If we cannot
        find a ring to sort an atom into, we will create a new ring containing that atom.
        It is important that we always check atoms that are already connected to existing rings
        before completely disconnected atoms. Otherwise, we will erroneously create new rings.
        """
        cbAtomList = copyGroup.sortByConnectivity(cbAtomList)
        cbfAtomList1 = copyGroup.sortByConnectivity(cbfAtomList1)
        cbfAtomList2 = copyGroup.sortByConnectivity(cbfAtomList2)

        """
        Step 4. Initalize the list of rings with any benzene rings that are already explicitly stated

        The variable rings is a list of lists. Each list in rings represents one full benzene rings,
        so it will eventually have six benzene carbons in it. Each ring's list will have the atoms
        sorted by connectivity, such that any atom is bonded to the atoms preceding and following it
        in the list. The first and last atom of the list will also be bonded together.

        """
        rings=[cycle for cycle in copyGroup.getAllCyclesOfSize(6) if Group(atoms = cycle).isAromaticRing()]

        """
        Step 5. Add Cbf2 atoms to the correct rings

        In this step, we define 'ring seeds', three atom combination unique to a each benzne ring. We
        then try to merge each ring seed into existing rings. If nosuittable ring is found, we create
        a new ring from the ring seed.

        Every Cbf2 atom is in two rings with unique ring seeds defined by partneredCbf-Cbf2atom-otherCbf
        and partneredCbf-Cbf2Atom-CbAtom. We may need to create the Cbatom in the last ring seed if it
        is not available.
        """
        for cbfAtom in cbfAtomList2:
            if connectedCbfs[cbfAtom][0] is partners[cbfAtom]: otherCbf = connectedCbfs[cbfAtom][1]
            else: otherCbf = connectedCbfs[cbfAtom][0]
            #These two ring seeds represent the two unique rings
            newRingSeeds = [[partners[cbfAtom], cbfAtom, otherCbf],
                            [partners[cbfAtom], cbfAtom]]
            allLigands = cbfAtom.bonds.keys()
            #add a new cb atom to the second newRing seed
            if len(allLigands) == 2:
                newAtom = copyGroup.createAndConnectAtom(['Cb'], cbfAtom, [1.5])
                newRingSeeds[1].append(newAtom)
            #join the existing atom to the ringSeed
            elif len(allLigands) == 3:
                for atom2 in allLigands:
                    if atom2 not in connectedCbfs[cbfAtom]:
                        newRingSeeds[1].append(atom2)
                        break
            #Check for duplicates, merge or create new rings
            for index1, ring1 in enumerate(newRingSeeds):
                mergeRingDict={}
                for index, ring2 in enumerate(rings):
                    #check if a duplicate of a fully created ring
                    if checkSet(ring2, ring1): break
                    #Next try to merge the ringseed into rings
                    mergeRing = mergeOverlappingBenzeneRings(ring2, ring1, 2)
                    if mergeRing:
                        mergeRingDict[index] = mergeRing
                        break
                #otherwise add this ringSeed because it represents a completely new ring
                else: rings.append(ring1)
                #if we merged a ring, we need to remove the old ring from rings and add the merged ring
                rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]

        """
        Step 6. Add Cbf1 atoms to the correct rings

        Every Cbf1 atom is in two rings with its partner. In this step, we add this ring seed
        twice to the rings.
        """
        for cbfAtom in cbfAtomList1:
            newRingSeed = [partners[cbfAtom], cbfAtom]
            inRing = 0
            #check to see if duplicate of an existing ring
            for ring in rings:
                if checkSet(ring, newRingSeed):
                    inRing +=1
            #move on to next cbfAtom if we found two rings
            if inRing ==2: continue
            #try to merge into existing rings, if cbf1 is connected
            for index, ring in enumerate(rings):
                mergeRingDict={}
                mergeRing = mergeOverlappingBenzeneRings(ring, newRingSeed, 1)
                if mergeRing:
                    inRing+=1
                    mergeRingDict[index]=mergeRing
                    #move on to next cbfAtom if we found two rings
                    if inRing ==2: break
            #if we merged a ring, we need to remove the old ring from rings and add the merged ring
            rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]
            #if we still dont have two ring, we create a completely new ring
            if inRing < 2:
                for x in range(2-inRing):
                    rings.append(copy(newRingSeed))

        """
        Step 7. Add Cb atoms to the correct rings

        Each Cb atom is part of exactly one benzene ring. In this step we merge or make
        a new ring for each Cb atom.
        """
        for cbAtom in cbAtomList:
            inRing = 0
            #check to see if already in a ring
            for ring in rings:
                if checkSet(ring, [cbAtom]): inRing +=1
            #move on to next ring cbAtom if in a ring
            if inRing == 1 : continue
            #check to see if can be merged to an existing ring
            for index, ring in enumerate(rings):
                mergeRingDict={}
                mergeRing = addCbAtomToRing(ring, cbAtom)
                if mergeRing:
                    inRing+=1
                    mergeRingDict[index]=mergeRing
                    break
            #if we merged a ring, we need to remove the old ring from rings and add the merged ring
            rings = [rings[index] if not index in mergeRingDict else mergeRingDict[index] for index in range(len(rings))]
            #Start completely new ring if not any of above true
            if inRing == 0:
                rings.append([cbAtom])

        """
        Step 8. Grow each partial ring up to six carbon atoms

        In this step we create new Cb atoms and add them to any rings which do not have 6 atoms.
        """
        mergedRingDict={}
        for index, ring in enumerate(rings):
            carbonsToGrow = 6-len(ring)
            mergedRingDict[index] = []
            for x in range(carbonsToGrow):
                if x ==0: lastAtom = ring[-1]
                else: lastAtom = mergedRingDict[index][-1]
                #add a new atom to the ring and the group
                newAtom = copyGroup.createAndConnectAtom(['Cb'], lastAtom, [1.5])
                mergedRingDict[index].append(newAtom)
                #At the end attach to the other endpoint
                if x == carbonsToGrow -1:
                    newBond = GroupBond(ring[0], newAtom, order=[1.5])
                    copyGroup.addBond(newBond)

        return copyGroup

    def pickWildcards(self):
        """
        Returns: the :class:Group object without wildcards in either atomtype or bonding

        This function will naively pick the first atomtype for each atom, but will try
        to pick bond orders that make sense given the selected atomtypes
        """
        for atom1 in self.atoms:
            atom1.atomType=[atom1.atomType[0]]
            for atom2, bond12 in atom1.bonds.iteritems():
                #skip dynamic bond ordering if there are no wildcards
                if len(bond12.order) < 2 :
                    continue
                atom1Features = atom1.atomType[0].getFeatures()
                atom2Features = atom2.atomType[0].getFeatures()
                #skip dynamic bond ordering if there are no features required by the atomtype
                if not any(atom1Features) and not any(atom2Features):
                    bond12.order = [bond12.order[0]]
                    atom2.bonds[atom1].order = bond12.order
                    continue

                atom1Bonds = atom1.countBonds()
                atom2Bonds = atom2.countBonds()
                requiredFeatures1 = [atom1Features[x][0] - atom1Bonds[x] if atom1Features[x] else 0 for x in range(len(atom1Bonds))]
                requiredFeatures2 = [atom2Features[x][0] - atom2Bonds[x] if atom2Features[x] else 0 for x in range(len(atom2Bonds))]

                #subtract 1 from allDouble for each sDouble, oDouble, rDouble so that we don't count them twice
                requiredFeatures1[1] = requiredFeatures1[1] - requiredFeatures1[2] - requiredFeatures1[3] -requiredFeatures1[4]
                requiredFeatures2[1] = requiredFeatures2[1] - requiredFeatures2[2] - requiredFeatures2[3] -requiredFeatures2[4]
                #reverse it because coincidentally the reverse has good priority on what to check first
                requiredFeatures1.reverse()
                requiredFeatures2.reverse()

                #required features are a now list of [benzene, triple, sDouble, oDouble, rDouble, allDouble, single]
                for index, (feature1, feature2) in enumerate(zip(requiredFeatures1[:-1], requiredFeatures2[:-1])):
                    if feature1 > 0 or feature2 > 0:
                        if index == 0 and 1.5 in bond12.order: #benzene bonds
                            bond12.order = [1.5]
                            atom2.bonds[atom1].order = bond12.order
                            break
                        elif index == 1 and 3 in bond12.order: #triple bond
                            bond12.order = [3]
                            atom2.bonds[atom1].order = bond12.order
                            break
                        elif index > 1 and 2 in bond12.order: #any case of double bonds
                            if index == 2: #sDouble bonds
                                if (feature1 > 0 and atom2.isSulfur()) or (feature2 > 0 and atom1.isSulfur()):
                                    bond12.order = [2]
                                    atom2.bonds[atom1].order = bond12.order
                                    break
                            elif index == 3: #oDoubleBonds
                                if (feature1 > 0 and atom2.isOxygen()) or (feature2 > 0 and atom1.isOxygen()):
                                    bond12.order = [2]
                                    atom2.bonds[atom1].order = bond12.order
                                    break
                            else: #rDouble or allDouble necessary
                                bond12.order = [2]
                                atom2.bonds[atom1].order = bond12.order
                                break
                else: #no features required, then pick the first order
                    bond12.order = [bond12.order[0]]
                    atom2.bonds[atom1].order = bond12.order

        #if we have wildcard atomtypes pick one based on ordering of allElements
        for atom in self.atoms:
            for elementLabel in allElements:
                if atomTypes[elementLabel] in atom.atomType[0].specific:
                    atom.atomType=[atomTypes[elementLabel]]
                    break

    def makeSampleMolecule(self):
        """
        Returns: A sample class :Molecule: from the group
        """

        modifiedGroup = self.copy(deep = True)

        #Remove all wildcards
        modifiedGroup.pickWildcards()

        #check that there are less than three Cbf atoms
        cbfCount = 0
        for atom in modifiedGroup.atoms:
            if atom.atomType[0] is atomTypes['Cbf']: cbfCount+=1
        if cbfCount > 3:
            if not modifiedGroup.isBenzeneExplicit():
                raise ImplicitBenzeneError("{0} has more than three Cbf atoms and does not have fully explicit benzene rings.")

        #Add implicit atoms
        modifiedGroup = modifiedGroup.addImplicitAtomsFromAtomType()

        #Add implicit benzene rings
        if not modifiedGroup.isBenzeneExplicit():
            modifiedGroup = modifiedGroup.addImplicitBenzene()
        #Make dictionary of :GroupAtoms: to :Atoms: and vice versa
        groupToMol = {}
        molToGroup = {}
        for atom in modifiedGroup.atoms:
            molAtom = atom.makeSampleAtom()
            groupToMol[atom] = molAtom
            molToGroup[molAtom] = atom

        #create the molecule
        newMolecule = mol.Molecule(atoms = groupToMol.values())
            
        #Add explicit bonds to :Atoms:
        for atom1 in modifiedGroup.atoms:
            for atom2, bond12 in atom1.bonds.iteritems():
                bond12.makeBond(newMolecule, groupToMol[atom1], groupToMol[atom2])

        #Saturate up to expected valency
        for molAtom in newMolecule.atoms:
            if molAtom.charge:
                statedCharge = molAtom.charge
            #otherwise assume no charge (or implicit atoms we assume hvae no charge)
            else:
                statedCharge = 0
            molAtom.updateCharge()
            if molAtom.charge - statedCharge:
                hydrogenNeeded = molAtom.charge - statedCharge
                if molAtom in molToGroup and molToGroup[molAtom].atomType[0].single:
                    maxSingle = max(molToGroup[molAtom].atomType[0].single)
                    singlePresent = sum([1 for atom in molAtom.bonds if molAtom.bonds[atom].isSingle()])
                    maxHydrogen = maxSingle - singlePresent
                    if hydrogenNeeded > maxHydrogen: hydrogenNeeded = maxHydrogen
                for x in range(hydrogenNeeded):
                    newH = mol.Atom('H', radicalElectrons=0, lonePairs=0, charge=0)
                    newBond = mol.Bond(molAtom, newH, 1)
                    newMolecule.addAtom(newH)
                    newMolecule.addBond(newBond)
                molAtom.updateCharge()

        newMolecule.update()

        #Check that the charge of atoms is expected
        for atom in newMolecule.atoms:
            if abs(atom.charge) > 0:
                if atom in molToGroup:
                    groupAtom = molToGroup[atom]
                else:
                    raise UnexpectedChargeError(graph = newMolecule)
                #check hardcoded atomtypes
                positiveCharged = ['Csc','Cdc',
                                   'N3sc','N5sc','N5dc','N5ddc','N5tc','N5b',
                                   'O2sc','O4sc','O4dc','O4tc',
                                   'S2sc','S4sc','S4dc','S4tdc','S6sc','S6dc','S6tdc']
                negativeCharged = ['C2sc','C2dc','C2tc',
                                   'N0sc','N1sc','N1dc','N5dddc',
                                   'O0sc',
                                   'S0sc','S2sc','S2dc','S2tc','S4dc','S4tdc','S6sc','S6dc','S6tdc']
                if groupAtom.atomType[0] in [atomTypes[x] for x in positiveCharged] and atom.charge > 0:
                    pass
                elif groupAtom.atomType[0] in [atomTypes[x] for x in negativeCharged] and atom.charge < 0:
                    pass
                #declared charge in original group is not same as new charge
                elif atom.charge in groupAtom.charge:
                    pass
                else:
                    raise UnexpectedChargeError(graph = newMolecule)

        return newMolecule

    def isBenzeneExplicit(self):
        """

        Returns: 'True' if all Cb, Cbf atoms are in completely explicitly stated benzene rings.

        Otherwise return 'False'

        """

        #classify atoms
        cbAtomList = []

        #only want to work with carbon atoms
        labelsOfCarbonAtomTypes = [x.label for x in atomTypes['C'].specific] + ['C', 'N3b', 'N5b']

        for atom in self.atoms:
            if not atom.atomType[0].label in labelsOfCarbonAtomTypes: continue
            elif atom.atomType[0].label in ['Cb', 'Cbf', 'N3b', 'N5b']: #Make Cb and N3b into normal cb atoms
                cbAtomList.append(atom)
            else:
                benzeneBonds = 0
                for atom2, bond12 in atom.bonds.iteritems():
                    if bond12.isBenzene(): benzeneBonds+=1
                if benzeneBonds >0: cbAtomList.append(atom)

        #get all explicit benzene rings
        rings=[cycle for cycle in self.getAllCyclesOfSize(6) if Group(atoms = cycle).isAromaticRing()]

        #test that all benzene atoms are in benzene rings
        for atom in cbAtomList:
            inRing = False
            for ring in rings:
                if atom in ring: inRing = True
            if not inRing: return False
        else: return True

    def mergeGroups(self, other):
        """
        This function takes `other` :class:Group object and returns a merged :class:Group object based
        on overlapping labeled atoms between self and other

        Currently assumes `other` can be merged at the closest labelled atom
        """
        labeled1 = self.getLabeledAtoms()
        labeled2 = other.getLabeledAtoms()
        overlappingLabels = [x for x in labeled1 if x in labeled2]

        #dictionary of key = original atoms, value = copy atoms for deep copies of the two groups
        selfDict= self.copyAndMap()
        otherDict = other.copyAndMap()

        #Sort atoms to go into the mergedGroup
        mergedGroupAtoms = otherDict.values() #all end atoms with end up in the mergedGroup
        #only non-overlapping atoms from the backbone will be in the mergedGroup
        for originalAtom, newAtom in selfDict.iteritems():
            if not originalAtom.label in overlappingLabels:
                mergedGroupAtoms.append(newAtom)

        mergedGroup = Group(atoms=mergedGroupAtoms)

        """
        The following loop will move bonds that are exclusively in the new backbone so that
        they connect to the backbone. For example, assume backbone has bond atomA-atomB,
        where atomB is labelled as *2 and there is an atomC in the end analgously labelled
        *2. We need to remove the bond between atomA and atomB. Then we need to add a bond
        between atomA and atomC.
        """
        bondsToRemove = []
        for label in overlappingLabels:
            oldAtomB = self.getLabeledAtom(label)
            for oldAtomA, oldBondAB in oldAtomB.bonds.iteritems():
                if not oldAtomA.label in overlappingLabels: #this is bond we need to transfer over
                    #find and record bondAB from new backbone for later removal
                    newAtomA = selfDict[oldAtomA]
                    newAtomB = selfDict[oldAtomB]
                    newAtomC = mergedGroup.getLabeledAtom(oldAtomB.label)
                    for atom, newBondAB in newAtomA.bonds.iteritems():
                        if atom is newAtomB:
                            bondsToRemove.append(newBondAB)
                            break
                    #add bond between atomA and AtomC
                    newBondAC = GroupBond(newAtomA, newAtomC, order= oldBondAB.order)
                    mergedGroup.addBond(newBondAC)
        #remove bonds from mergedGroup
        for bond in bondsToRemove:
            mergedGroup.removeBond(bond)

        return mergedGroup

    def resetRingMembership(self):
        """
        Resets ring membership information in the GroupAtom.props attribute.
        """
        cython.declare(ratom=GroupAtom)

        for atom in self.atoms:
            if 'inRing' in atom.props:
                del atom.props['inRing']
