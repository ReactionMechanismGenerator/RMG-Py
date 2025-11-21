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
from rmgpy.molecule.atomtype import ATOMTYPES, allElements, nonSpecifics, get_features, AtomType
from rmgpy.molecule.element import PeriodicSystem
from rmgpy.molecule.graph import Vertex, Edge, Graph
from rmgpy.molecule.fragment import CuttingLabel

# helper functions
# these were originally nested inside the indicated parent function, but when we upgraded to
# Cython 3 this was no longer allowed - thus, they now live here.

# add_implicit_benzene
def check_set(super_list, sub_list):
    """
    Args:
        super_list: list to check if superset of partList
        sub_list:  list to check if subset of superList

    Returns: Boolean to see if super_list is a superset of sub_list

    """
    super_set = set(super_list)
    sub_set = set(sub_list)
    return super_set.issuperset(sub_set)

def add_cb_atom_to_ring(ring, cb_atom):
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
        if bond12.is_benzene():
            if atom2 is ring[-1]:
                merged_ring = ring + [cb_atom]
            elif atom2 is ring[0]:
                merged_ring = [cb_atom] + ring

    return merged_ring

def merge_overlapping_benzene_rings(ring1, ring2, od):
    """
    The input arguments of rings are always in the order that the atoms appear
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

class GroupAtom(Vertex):
    """
    An atom group. This class is based on the :class:`Atom` class, except that
    it uses :ref:`atom types <atom-types>` instead of elements, and all
    attributes are lists rather than individual values. The attributes are:

    ==================== =================== ====================================
    Attribute            Type                Description
    ==================== =================== ====================================
    `atomtype`           ``list``            The allowed atom types (as :class:`AtomType` objects)
    `radical_electrons`  ``list``            The allowed numbers of radical electrons (as short integers)
    `charge`             ``list``            The allowed formal charges (as short integers)
    `label`              ``str``             A string label that can be used to tag individual atoms
    `lone_pairs`         ``list``            The number of lone electron pairs
    `charge`             ``list``            The partial charge of the atom
    `site`               ``list``            The allowed adsorption sites
    `morphology`         ``list``            The allowed morphologies
    `props`              ``dict``            Dictionary for storing additional atom properties
    `reg_dim_atm`        ``list``            List of atom types that are free dimensions in tree optimization
    `reg_dim_u`          ``list``            List of unpaired electron numbers that are free dimensions in tree optimization
    `reg_dim_r`          ``list``            List of inRing values that are free dimensions in tree optimization
    `reg_dim_site`       ``list``            List of sites that are free dimensions in tree optimization
    `reg_dim_morphology` ``list``            List of morphologies that are free dimensions in tree optimization
    ==================== =================== ====================================

    Each list represents a logical OR construct, i.e. an atom will match the
    group if it matches *any* item in the list. However, the
    `radical_electrons`, and `charge` attributes are linked
    such that an atom must match values from the same index in each of these in
    order to match.
    """

    def __init__(self, atomtype=None, radical_electrons=None, charge=None, label='', lone_pairs=None, site=None, morphology=None,
                 props=None):
        Vertex.__init__(self)
        self.atomtype = atomtype or []
        for index in range(len(self.atomtype)):
            if isinstance(self.atomtype[index], str):
                self.atomtype[index] = ATOMTYPES[self.atomtype[index]]
        self.radical_electrons = radical_electrons or []
        self.charge = charge or []
        self.label = label
        self.lone_pairs = lone_pairs or []
        self.site = site or []
        self.morphology = morphology or []
        self.props = props or {}

        self.reg_dim_atm = [[], []]
        self.reg_dim_u = [[], []]
        self.reg_dim_r = [[], []]
        self.reg_dim_site = [[], []]
        self.reg_dim_morphology = [[], []]

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        d = {
            'edges': self.edges,
            'connectivity1': self.connectivity1,
            'connectivity2': self.connectivity2,
            'connectivity3': self.connectivity3,
            'sorting_label': self.sorting_label,
        }
        atomtype = self.atomtype
        if atomtype is not None:
            atomtype = [a.label for a in atomtype]
        return (GroupAtom, (atomtype, self.radical_electrons, self.charge, self.label, self.lone_pairs, self.site,
                            self.morphology, self.props), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sorting_label = d['sorting_label']

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '[{0} {1}]'.format(self.label, ','.join([repr(a.label) for a in self.atomtype]))

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
            self.atomtype[:],
            self.radical_electrons[:],
            self.charge[:],
            self.label,
            self.lone_pairs[:],
            self.site[:],
            self.morphology[:],
            deepcopy(self.props),
        )

    def _change_bond(self, order):
        """
        Update the atom group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        atomtype = []
        for atom in self.atomtype:
            if order == 1:
                atomtype.extend(atom.increment_bond)
            elif order == -1:
                atomtype.extend(atom.decrement_bond)
            else:
                raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: '
                                  'Invalid order "{0}".'.format(order))
        if len(atomtype) == 0:
            raise ActionError('Unable to update GroupAtom due to CHANGE_BOND action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))
        # Set the new atom types, removing any duplicates
        self.atomtype = list(set(atomtype))

    def _form_bond(self, order):
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
        atomtype = []
        for atom in self.atomtype:
            atomtype.extend(atom.form_bond)
        if len(atomtype) == 0:
            raise ActionError('Unable to update GroupAtom due to FORM_BOND action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))
        # Set the new atom types, removing any duplicates
        self.atomtype = list(set(atomtype))

    def _break_bond(self, order):
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
        atomtype = []
        for atom in self.atomtype:
            atomtype.extend(atom.break_bond)
        if len(atomtype) == 0:
            raise ActionError('Unable to update GroupAtom due to BREAK_BOND action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))
        # Set the new atom types, removing any duplicates
        self.atomtype = list(set(atomtype))

    def _gain_radical(self, radical):
        """
        Update the atom group as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.

        The 'radicalElectron' attribute can be an empty list if we use the wildcard
        argument ux in the group definition. In this case, we will have this
        function set the atom's 'radicalElectron' to a list allowing from `radical`
        up to 4 radical electrons.
        """
        radical_electrons = []
        if any([len(atomtype.increment_radical) == 0 for atomtype in self.atomtype]):
            raise ActionError('Unable to update GroupAtom due to GAIN_RADICAL action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))
        if not self.radical_electrons:
            radical_electrons = [1, 2, 3, 4][radical-1:]
        else:
            for electron in self.radical_electrons:
                radical_electrons.append(electron + radical)
        # Set the new radical electron counts
        self.radical_electrons = radical_electrons

    def _lose_radical(self, radical):
        """
        Update the atom group as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.

        The 'radicalElectron' attribute can be an empty list if we use the wildcard
        argument ux in the group definition. In this case, we will have this
        function set the atom's 'radicalElectron' to a list allowing 0, 1, 2,
        or 3 radical electrons.
        """
        radical_electrons = []
        if any([len(atomtype.decrement_radical) == 0 for atomtype in self.atomtype]):
            raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))

        if not self.radical_electrons:
            radical_electrons = [0, 1, 2, 3]
        else:
            for electron in self.radical_electrons:
                electron = electron - radical
                if electron < 0:
                    raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: '
                                      'Invalid radical electron set "{0}".'.format(self.radical_electrons))
                radical_electrons.append(electron)

        # Set the new radical electron counts
        self.radical_electrons = radical_electrons

    def _gain_charge(self, charge):
        """
        Update the atom group as a result of applying a GAIN_CHARGE action,
        where `charge` specifies the charge gained.
        """
        atomtype = []

        for atom in self.atomtype:
            atomtype.extend(atom.increment_charge)

        if any([len(atom.increment_charge) == 0 for atom in self.atomtype]):
            raise ActionError('Unable to update GroupAtom due to GAIN_CHARGE action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))

        if isinstance(self.charge,list):
            charges = []
            for c in self.charge:
                charges.append(c+charge)
            self.charge = charges
        else:
            self.charge += 1

        self.atomtype = list(set(atomtype))

    def _lose_charge(self, charge):
        """
        Update the atom group as a result of applying a LOSE_CHARGE action,
        where `charge` specifies lost charge.
        """
        atomtype = []

        for atom in self.atomtype:
            atomtype.extend(atom.decrement_charge)

        if any([len(atomtype.decrement_charge) == 0 for atomtype in self.atomtype]):
            raise ActionError('Unable to update GroupAtom due to LOSE_CHARGE action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))

        if isinstance(self.charge,list):
            charges = []
            for c in self.charge:
                charges.append(c-charge)
            self.charge = charges
        else:
            self.charge -= 1

        self.atomtype = list(set(atomtype))

    def _gain_pair(self, pair):
        """
        Update the atom group as a result of applying a GAIN_PAIR action,
        where `pair` specifies the number of lone electron pairs to add.
        """
        lone_pairs = []
        atomtype = []

        for atom in self.atomtype:
            atomtype.extend(atom.increment_lone_pair)
        if any([len(atom.increment_lone_pair) == 0 for atom in self.atomtype]):
            raise ActionError('Unable to update GroupAtom due to GAIN_PAIR action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))

        # Add a lone pair to a group atom with none
        if not self.lone_pairs:
            self.lone_pairs = [1, 2, 3, 4][pair-1:]  # set to a wildcard of any number greater than or equal to `pair`
        # Add a lone pair to a group atom that already has at least one lone pair
        else:
            for x in self.lone_pairs:
                lone_pairs.append(x + pair)
            # Set the new lone electron pair count
            self.lone_pairs = lone_pairs

        # Set the new atom types, removing any duplicates
        self.atomtype = list(set(atomtype))

    def _lose_pair(self, pair):
        """
        Update the atom group as a result of applying a LOSE_PAIR action,
        where `pair` specifies the number of lone electron pairs to remove.
        """
        lone_pairs = []
        atomtype = []

        for atom in self.atomtype:
            atomtype.extend(atom.decrement_lone_pair)
        if any([len(atom.decrement_lone_pair) == 0 for atom in self.atomtype]):
            raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: '
                              'Unknown atom type produced from set "{0}".'.format(self.atomtype))

        if not self.lone_pairs:
            self.lone_pairs = [0, 1, 2, 3]  # set to a wildcard of any number fewer than 4
        else:
            for x in self.lone_pairs:
                if x - pair < 0:
                    raise ActionError('Unable to update GroupAtom due to LOSE_PAIR action: '
                                      'Invalid lone electron pairs set "{0}".'.format(self.lone_pairs))
                lone_pairs.append(x - pair)
            # Set the new lone electron pair count
            self.lone_pairs = lone_pairs

        # Set the new atom types, removing any duplicates
        self.atomtype = list(set(atomtype))

    def apply_action(self, action):
        """
        Update the atom group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        act = action[0].upper()
        if act == 'CHANGE_BOND':
            self._change_bond(action[2])
        elif act == 'FORM_BOND':
            self._form_bond(action[2])
        elif act == 'BREAK_BOND':
            self._break_bond(action[2])
        elif act == 'GAIN_RADICAL':
            self._gain_radical(action[2])
        elif act == 'GAIN_CHARGE':
            self._gain_charge(action[2])
        elif act == 'LOSE_RADICAL':
            self._lose_radical(action[2])
        elif act == 'LOSE_CHARGE':
            self._lose_charge(action[2])
        elif action[0].upper() == 'GAIN_PAIR':
            self._gain_pair(action[2])
        elif action[0].upper() == 'LOSE_PAIR':
            self._lose_pair(action[2])
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
        for atomType1 in self.atomtype:
            for atomType2 in group.atomtype:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        for atomType1 in group.atomtype:
            for atomType2 in self.atomtype:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        for radical1 in self.radical_electrons:
            if group.radical_electrons:  # Only check if the list is non-empty.  An empty list indicates a wildcard.
                for radical2 in group.radical_electrons:
                    if radical1 == radical2: break
                else:
                    return False
        for radical1 in group.radical_electrons:
            if self.radical_electrons:
                for radical2 in self.radical_electrons:
                    if radical1 == radical2: break
                else:
                    return False
        for lp1 in self.lone_pairs:
            if group.lone_pairs:
                for lp2 in group.lone_pairs:
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
        # Each site in self must have an equivalent in other (and vice versa)
        for site1 in self.site:
            if group.site:
                for site2 in group.site:
                    if site1 == site2: break
                else:
                    return False
        for site1 in group.site:
            if self.site:
                for site2 in self.site:
                    if site1 == site2: break
                else:
                    return False
        # Each morphology in self must have an equivalent in other (and vice versa)
        for morphology1 in self.morphology:
            if group.morphology:
                for morphology2 in group.morphology:
                    if morphology1 == morphology2: break
                else:
                    return False
        for morphology1 in group.morphology:
            if self.morphology:
                for morphology2 in self.morphology:
                    if morphology1 == morphology2: break
                else:
                    return False
        # Other properties must have an equivalent in other (and vice versa)
        # Absence of the 'inRing' prop indicates a wildcard
        if 'inRing' in self.props and 'inRing' in group.props:
            if self.props['inRing'] != group.props['inRing']:
                return False
        # Otherwise the two atom groups are equivalent
        return True

    def is_specific_case_of(self, other):
        """
        Returns ``True`` if `self` is the same as `other` or is a more
        specific case of `other`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """
        cython.declare(group=GroupAtom)
        if not isinstance(other, GroupAtom):
            # Let the is_specific_case_of method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.is_specific_case_of(self)
        group = other

        cython.declare(atomType1=AtomType, atomtype2=AtomType, radical1=cython.short, radical2=cython.short,
                       lp1=cython.short, lp2=cython.short, charge1=cython.short, charge2=cython.short,
                       site1=str, site2=str, morphology1=str, morphology2=str)
        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomtype:  # all these must match
            for atomType2 in group.atomtype:  # can match any of these
                if atomType1.is_specific_case_of(atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        if self.radical_electrons:
            for radical1 in self.radical_electrons:
                if group.radical_electrons:
                    for radical2 in group.radical_electrons:
                        if radical1 == radical2: break
                    else:
                        return False
        else:
            if group.radical_electrons: return False
        if self.lone_pairs:
            for lp1 in self.lone_pairs:
                if group.lone_pairs:
                    for lp2 in group.lone_pairs:
                        if lp1 == lp2: break
                    else:
                        return False
        else:
            if group.lone_pairs: return False
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
        # Each site in self must have an equivalent in other
        if self.site:
            for site1 in self.site:
                if group.site:
                    for site2 in group.site:
                        if site1 == site2: break
                    else:
                        return False
        else:
            if group.site: return False
        # Each morphology in self must have an equivalent in other
        if self.morphology:
            for morphology1 in self.morphology:
                if group.morphology:
                    for morphology2 in group.morphology:
                        if morphology1 == morphology2: break
                    else:
                        return False
        else:
            if group.morphology: return False
        # Other properties must have an equivalent in other
        # Absence of the 'inRing' prop indicates a wildcard
        if 'inRing' in self.props and 'inRing' in group.props:
            if self.props['inRing'] != group.props['inRing']:
                return False
        elif 'inRing' not in self.props and 'inRing' in group.props:
            return False
        # Otherwise self is in fact a specific case of other
        return True

    def is_surface_site(self):
        """
        Return ``True`` if the atom represents a surface site or ``False`` if not.
        """
        site_type = ATOMTYPES['X']
        return all([s.is_specific_case_of(site_type) for s in self.atomtype])

    def is_bonded_to_surface(self):
        """
        Return ``True`` if the atom is bonded to a surface GroupAtom `X`
        ``False`` if it is not
        """
        cython.declare(bonded_atom=GroupAtom)
        for bonded_atom in self.bonds.keys():
            if bonded_atom.is_surface_site():
                return True
        return False

    def is_electron(self):
        """
        Return ``True`` if the atom represents a surface site or ``False`` if not.
        """
        return self.atomtype[0] == ATOMTYPES['e']

    def is_proton(self):
        """
        Return ``True`` if the atom represents a surface site or ``False`` if not.
        """
        return self.atomtype[0] == ATOMTYPES['H+']

    def is_oxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if not.
        """
        all_oxygens = [ATOMTYPES['O']] + ATOMTYPES['O'].specific
        check_list = [x in all_oxygens for x in self.atomtype]
        return all(check_list)

    def is_sulfur(self):
        """
        Return ``True`` if the atom represents an sulfur atom or ``False`` if not.
        """
        all_sulfur = [ATOMTYPES['S']] + ATOMTYPES['S'].specific
        check_list = [x in all_sulfur for x in self.atomtype]
        return all(check_list)

    def is_nitrogen(self):
        """
        Return ``True`` if the atom represents a nitrogen atom or ``False`` if not.
        """
        all_nitrogen = [ATOMTYPES['N']] + ATOMTYPES['N'].specific
        check_list = [x in all_nitrogen for x in self.atomtype]
        return all(check_list)

    def is_carbon(self):
        """
        Return ``True`` if the atom represents a carbon atom or ``False`` if not.
        """
        all_carbon = [ATOMTYPES['C']] + ATOMTYPES['C'].specific
        check_list = [x in all_carbon for x in self.atomtype]
        return all(check_list)

    def is_fluorine(self):
        """
        Return ``True`` if the atom represents a fluorine atom or ``False`` if not.
        """
        all_fluorine = [ATOMTYPES['F']] + ATOMTYPES['F'].specific
        check_list = [x in all_fluorine for x in self.atomtype]
        return all(check_list)

    def is_chlorine(self):
        """
        Return ``True`` if the atom represents a chlorine atom or ``False`` if not.
        """
        all_chlorine = [ATOMTYPES['Cl']] + ATOMTYPES['Cl'].specific
        check_list = [x in all_chlorine for x in self.atomtype]
        return all(check_list)

    def is_bromine(self):
        """
        Return ``True`` if the atom represents a bromine atom or ``False`` if not.
        """
        all_bromine = [ATOMTYPES['Br']] + ATOMTYPES['Br'].specific
        check_list = [x in all_bromine for x in self.atomtype]
        return all(check_list)

    def is_lithium(self):
        """
        Return ``True`` if the atom represents a bromine atom or ``False`` if not.
        """
        all_lithium = [ATOMTYPES['Li']] + ATOMTYPES['Li'].specific
        check_list = [x in all_lithium for x in self.atomtype]
        return all(check_list)

    def has_wildcards(self):
        """
        Return ``True`` if the atom has wildcards in any of the attributes:
        atomtype, radical electrons, lone pairs, charge, and bond order. Returns
        ''False'' if no attribute has wildcards.
        """
        if len(self.atomtype) > 1:
            return True
        elif len(self.radical_electrons) > 1 or len(self.radical_electrons) == 0:
            return True
        elif len(self.lone_pairs) > 1:
            return True
        for bond in self.bonds.values():
            if len(bond.order) > 1:
                return True
        return False

    def count_bonds(self, wildcards=False):
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
            if bond12.is_single(wildcards=True):
                single += 1
            if bond12.is_double(wildcards=True):
                if atom2.is_oxygen():
                    o_double += 1
                elif atom2.is_sulfur():
                    s_double += 1
                else:
                    # r_double is for double bonds NOT to oxygen or Sulfur
                    r_double += 1
            if bond12.is_triple(wildcards=True): triple += 1
            if bond12.is_quadruple(wildcards=True): quadruple += 1
            if bond12.is_benzene(wildcards=True): benzene += 1

        all_double = r_double + o_double + s_double

        # Warning: some parts of code assume this matches precisely the list returned by get_features()
        return [single, all_double, r_double, o_double, s_double, triple, quadruple, benzene]

    def make_sample_atom(self):
        """

        Returns: a class :Atom: object analagous to the GroupAtom

        This makes a sample, so it takes the first element when there are multiple options inside of
        self.atomtype, self.radical_electrons, self.lone_pairs, and self.charge

        """

        # Use the first atomtype to determine element, even if there is more than one atomtype
        atomtype = self.atomtype[0]
        element = None

        default_lone_pairs = {'H': 0,
                              'D': 0,
                              'T': 0,
                              'He': 1,
                              'C': 0,
                              'O': 2,
                              'N': 1,
                              'Si': 0,
                              'P': 1,
                              'S': 2,
                              'Ne': 4,
                              'Cl': 3,
                              'F': 3,
                              'Br': 3,
                              'I': 3,
                              'Ar': 4,
                              'X': 0,
                              'e': 0
                              }

        for element_label in allElements:
            if atomtype is ATOMTYPES[element_label] or atomtype in ATOMTYPES[element_label].specific:
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

        # Three possible values for charge and lone_pairs
        if self.charge:
            new_charge = self.charge[0]
        elif atomtype.charge:
            new_charge = atomtype.charge[0]
        else:
            new_charge = default_atom.charge

        if self.lone_pairs:
            new_lone_pairs = self.lone_pairs[0]
        elif atomtype.lone_pairs:
            new_lone_pairs = atomtype.lone_pairs[0]
        else:
            new_lone_pairs = default_atom.lone_pairs

        new_atom = mol.Atom(
            element=element,
            radical_electrons=self.radical_electrons[0] if self.radical_electrons else default_atom.radical_electrons,
            charge=new_charge,
            lone_pairs=new_lone_pairs,
            label=self.label if self.label else default_atom.label
        )

        # For some reason the default when no lone pairs is set to -100,
        # Based on git history, it is probably because RDKit requires a number instead of None
        if new_atom.lone_pairs == -100:
            new_atom.lone_pairs = default_lone_pairs[new_atom.symbol]

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
            self.set_order_str(order)
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

    def get_order_str(self):
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
            elif abs(value - 0.1) < 1e-4:
                values.append('H')
            elif abs(value - 0.05) < 1e-4:
                values.append('R')
            else:
                raise TypeError('Bond order number {} is not hardcoded as a string'.format(value))
        return values

    def set_order_str(self, new_order):
        """
        set the bond order using a valid bond-order character list
        """

        values = []
        for value in new_order:
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
            elif value == 'R':
                values.append(0.05)
            else:
                # try to see if an float disguised as a string was input by mistake
                try:
                    values.append(float(value))
                except ValueError:
                    raise TypeError('Bond order {} is not hardcoded into this method'.format(value))
        self.order = values

    def get_order_num(self):
        """
        returns the bond order as a list of numbers
        """
        return self.order

    def set_order_num(self, new_order):
        """
        change the bond order with a list of numbers
        """
        self.order = new_order

    def is_single(self, wildcards=False):
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

    def is_double(self, wildcards=False):
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

    def is_triple(self, wildcards=False):
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

    def is_quadruple(self, wildcards=False):
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

    def is_van_der_waals(self, wildcards=False):
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


    def is_reaction_bond(self, wildcards=False):
        """
        Return ``True`` if the bond represents a van der Waals bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are van der Waals.
        """
        if wildcards:
            for order in self.order:
                if abs(order[0] - 0.05) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0]-0.05) <= 1e-9 and len(self.order) == 1

    def is_benzene(self, wildcards=False):
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

    def is_hydrogen_bond(self, wildcards=False):
        """
        Return ``True`` if the bond represents a hydrogen bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are hydrogen bonds.
        """
        if wildcards:
            for order in self.order:
                if abs(order - 0.1) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 0.1) <= 1e-9 and len(self.order) == 1

    def is_reaction_bond(self, wildcards=False):
        """
        Return ``True`` if the bond represents a reaction bond or ``False`` if
        not. If `wildcards` is ``False`` we return False anytime there is more
        than one bond order, otherwise we return ``True`` if any of the options
        are reaction bonds.
        """
        if wildcards:
            for order in self.order:
                if abs(order - 0.05) <= 1e-9:
                    return True
            else:
                return False
        else:
            return abs(self.order[0] - 0.05) <= 1e-9 and len(self.order) == 1

    def _change_bond(self, order):
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

    def apply_action(self, action):
        """
        Update the bond group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self._change_bond(action[2])
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

    def is_specific_case_of(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """
        cython.declare(gb=GroupBond)
        if not isinstance(other, GroupBond):
            # Let the is_specific_case_of method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.is_specific_case_of(self)
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

    def make_bond(self, molecule, atom1, atom2):
        """
        Creates a :class: Bond between atom1 and atom2 analogous to self

        The intended input arguments should be class :Atom: not class :GroupAtom:
        Args:
            atom1: First :class: Atom the bond connects
            atom2: Second :class: Atom the bond connects

        """
        new_bond = mol.Bond(atom1, atom2, order=self.order[0])
        molecule.add_bond(new_bond)


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
    `metal`             ``list``            List of metals accepted for the group
    `facet`             ``list``            List of facets accepted for the group
    =================== =================== ====================================

    Corresponding alias methods to Molecule have also been provided.
    """

    def __init__(self, atoms=None, props=None, multiplicity=None, metal=None, facet=None):
        Graph.__init__(self, atoms)
        self.props = props or {}
        self.multiplicity = multiplicity or []
        self.metal = metal or []
        self.facet = facet or []
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

    def draw(self, file_format):
        """
        Use pydot to draw a basic graph of the group.

        Use format to specify the desired output file_format, eg. 'png', 'svg', 'ps', 'pdf', 'plain', etc.
        """
        import pydot

        graph = pydot.Dot(graph_type='graph', dpi="52")
        for index, atom in enumerate(self.atoms):
            atom_type = '{0!s} {1!s} '.format(index+1, atom.label if atom.label != '' else '')
            atom_type += ','.join([at.label for at in atom.atomtype])
            atom_type = '"' + atom_type + '"'
            graph.add_node(pydot.Node(name=str(index + 1), label=atom_type, fontname="Helvetica", fontsize="16"))
        for atom1 in self.atoms:
            for atom2, bond in atom1.bonds.items():
                index1 = self.atoms.index(atom1)
                index2 = self.atoms.index(atom2)
                if index1 < index2:
                    bond_type = ','.join([order for order in bond.get_order_str()])
                    bond_type = '"' + bond_type + '"'
                    graph.add_edge(pydot.Edge(src=str(index1 + 1), dst=str(index2 + 1),
                                              label=bond_type, fontname="Helvetica", fontsize="16"))

        img = graph.create(prog='neato', format=file_format)
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

    def add_atom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.add_vertex(atom)

    def add_bond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        return self.add_edge(bond)

    def get_bonds(self, atom):
        """
        Return a list of the bonds involving the specified `atom`.
        """
        return self.get_edges(atom)

    def get_bond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.get_edge(atom1, atom2)

    def has_atom(self, atom):
        """
        Returns ``True`` if `atom` is an atom in the graph, or ``False`` if
        not.
        """
        return self.has_vertex(atom)

    def has_bond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.has_edge(atom1, atom2)

    def contains_surface_site(self):
        """
        Returns ``True`` iff the group contains an 'X' surface site.
        """
        cython.declare(atom=GroupAtom)
        for atom in self.atoms:
            if atom.is_surface_site():
                return True
        return False

    def is_surface_site(self):
        """Returns ``True`` iff the group is nothing but a surface site 'X'."""
        return len(self.atoms) == 1 and self.atoms[0].is_surface_site()

    def get_surface_sites(self):
        """
        Get a list of surface site GroupAtoms in the group.
        Returns:
            List(GroupAtom): A list containing the surface site GroupAtoms in the molecule
        """
        cython.declare(atom=GroupAtom)
        return [atom for atom in self.atoms if atom.is_surface_site()]

    def is_proton(self):
        """Returns ``True`` iff the group is a proton"""
        return len(self.atoms) == 1 and self.atoms[0].is_proton()

    def is_electron(self):
        """Returns ``True`` iff the group is an electron"""
        return len(self.atoms) == 1 and self.atoms[0].is_electron()

    def remove_atom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        return self.remove_vertex(atom)

    def remove_bond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.remove_edge(bond)

    def remove_van_der_waals_bonds(self):
        """
        Remove all bonds that are definitely only van der Waals bonds.
        """
        cython.declare(bond=GroupBond)
        for bond in self.get_all_edges():
            if bond.is_van_der_waals(wildcards=False):
                self.remove_bond(bond)

    def sort_atoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        return self.sort_vertices()

    def sort_by_connectivity(self, atom_list):
        """
        Args:
            atom_list: input list of atoms

        Returns: a sorted list of atoms where each atom is connected to a previous
        atom in the list if possible
        """
        # if no input given just return
        if not atom_list: return atom_list

        sorted_atom_list = []
        sorted_atom_list.append(atom_list.pop(0))
        while atom_list:
            for atom1 in sorted_atom_list:
                added = False
                for atom2, bond12 in atom1.bonds.items():
                    if bond12.is_benzene() and atom2 in atom_list:
                        sorted_atom_list.append(atom2)
                        atom_list.remove(atom2)
                        added = True
                if added: break
            else:
                sorted_atom_list.append(atom_list.pop(0))

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

        self.update_connectivity_values()
        self.update_fingerprint()

    def update_charge(self):
        """
        Update the partial charge according to the valence electron, total bond order, lone pairs
        and radical electrons. This method is used for products of specific families with recipes that modify charges.
        """
        for atom in self.atoms:
            if isinstance(atom, CuttingLabel):
                continue
            if (len(atom.charge) == 1) and (len(atom.lone_pairs) == 1) and (len(atom.radical_electrons) == 1):
                # if the charge of the group is not labeled, then no charge update will be
                # performed. If there multiple charges are assigned, no update either.
                # Besides, this groupatom should have enough information to be updated
                atom_type = atom.atomtype[0]
                for element in allElements:
                    if atom_type is ATOMTYPES[element] or atom_type in ATOMTYPES[element].specific:
                        bond_order = 0
                        valence_electron = elements.PeriodicSystem.valence_electrons[element]
                        for _, bond in atom.bonds.items():
                            bond_order += bond.order[0]
                        lone_pairs = atom.lone_pairs[0]
                        radical_electrons = atom.radical_electrons[0]
                        atom.charge[0] = valence_electron - bond_order - 2 * lone_pairs - radical_electrons
                    else:
                        # if the group is not specified to specific element, charge will not be updated
                        pass

    def get_net_charge(self):
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

    def clear_reg_dims(self):
        """
        clear regularization dimensions
        """
        cython.declare(atm=GroupAtom)
        for atm in self.atoms:
            atm.reg_dim_atm = [[], []]
            atm.reg_dim_u = [[], []]
            atm.reg_dim_r = [[], []]
            atm.reg_dim_site = [[],[]]
            atm.reg_dim_morphology = [[],[]]
        for bd in self.get_all_edges():
            bd.reg_dim = [[], []]

    def get_extensions(self, r=None, r_bonds=None, r_un=None, basename='', atm_ind=None, atm_ind2=None, n_splits=None):
        """
        generate all allowed group extensions and their complements
        note all atomtypes except for elements and r/r!H's must be removed
        """
        cython.declare(atoms=list, atm=GroupAtom, atm2=GroupAtom, bd=GroupBond, i=int, j=int,
                       extents=list, RnH=list, typ=list)
        
        extents = []
        if r_bonds is None:
            r_bonds = [1, 1.5, 2, 3, 4]
        if r_un is None:
            r_un = [0, 1, 2, 3]
        if n_splits is None:
            n_splits = len(self.split())

        # generate appropriate r and r!H
        if r is None:
            r = elements.bde_elements  # set of possible r elements/atoms
            r = [ATOMTYPES[x] for x in r]

        R = r[:]
        if ATOMTYPES['X'] in R:
            R.remove(ATOMTYPES['X'])
        RnH = R[:]
        if ATOMTYPES['H'] in RnH:
            RnH.remove(ATOMTYPES['H'])

        atoms = self.atoms
        if atm_ind is None:
            for i, atm in enumerate(atoms):
                typ = atm.atomtype
                if atm.reg_dim_atm[0] == []:
                    if len(typ) == 1:
                        if typ[0].label == 'R':
                            extents.extend(self.specify_atom_extensions(i, basename, R))  # specify types of atoms
                        elif typ[0].label == 'R!H':
                            extents.extend(self.specify_atom_extensions(i, basename, RnH))
                        elif typ[0].label == 'Rx':
                            extents.extend(self.specify_atom_extensions(i, basename, r))

                    else:
                        extents.extend(self.specify_atom_extensions(i, basename, typ))
                else:
                    if len(typ) == 1:
                        if typ[0].label == 'R':
                            extents.extend(
                                self.specify_atom_extensions(i, basename, atm.reg_dim_atm[0]))  # specify types of atoms
                        elif typ[0].label == 'R!H':
                            extents.extend(
                                self.specify_atom_extensions(i, basename, list(set(atm.reg_dim_atm[0]) & set(R))))
                        elif typ[0].label == 'Rx':
                            extents.extend(
                                self.specify_atom_extensions(i, basename, list(set(atm.reg_dim_atm[0]) & set(r))))
                    else:
                        extents.extend(
                            self.specify_atom_extensions(i, basename, list(set(typ) & set(atm.reg_dim_atm[0]))))
                if atm.reg_dim_u[0] == []:
                    if len(atm.radical_electrons) != 1:
                        if len(atm.radical_electrons) == 0:
                            extents.extend(self.specify_unpaired_extensions(i, basename, r_un))
                        else:
                            extents.extend(self.specify_unpaired_extensions(i, basename, atm.radical_electrons))
                else:
                    if len(atm.radical_electrons) != 1 and len(atm.reg_dim_u[0]) != 1:
                        if len(atm.radical_electrons) == 0:
                            extents.extend(self.specify_unpaired_extensions(i, basename, atm.reg_dim_u[0]))
                        else:
                            extents.extend(self.specify_unpaired_extensions(i, basename, list(
                                set(atm.radical_electrons) & set(atm.reg_dim_u[0]))))
                if atm.reg_dim_r[0] == [] and 'inRing' not in atm.props:
                    extents.extend(self.specify_ring_extensions(i, basename))

                extents.extend(self.specify_external_new_bond_extensions(i, basename, r_bonds))
                for j, atm2 in enumerate(atoms):
                    if j < i and not self.has_bond(atm, atm2):
                        extents.extend(self.specify_internal_new_bond_extensions(i, j, n_splits, basename, r_bonds))
                    elif j < i:
                        bd = self.get_bond(atm, atm2)
                        if len(bd.order) > 1 and bd.reg_dim[0] == []:
                            extents.extend(self.specify_bond_extensions(i, j, basename, bd.order))
                        elif len(bd.order) > 1 and len(bd.reg_dim[0]) > 1 and len(bd.reg_dim[0]) > len(bd.reg_dim[1]):
                            extents.extend(self.specify_bond_extensions(i, j, basename, bd.reg_dim[0]))

        elif atm_ind is not None and atm_ind2 is not None:  # if both atm_ind and atm_ind2 are defined only look at the bonds between them
            i = atm_ind
            j = atm_ind2
            atm = atoms[i]
            atm2 = atoms[j]
            if j < i and not self.has_bond(atm, atm2):
                extents.extend(self.specify_internal_new_bond_extensions(i, j, n_splits, basename, r_bonds))
            if self.has_bond(atm, atm2):
                bd = self.get_bond(atm, atm2)
                if len(bd.order) > 1 and bd.reg_dim[0] == []:
                    extents.extend(self.specify_bond_extensions(i, j, basename, bd.order))
                elif len(bd.order) > 1 and len(bd.reg_dim[0]) > 1 and len(bd.reg_dim[0]) > len(bd.reg_dim[1]):
                    extents.extend(self.specify_bond_extensions(i, j, basename, bd.reg_dim[0]))

        elif atm_ind is not None:  # look at the atom at atm_ind
            i = atm_ind
            atm = atoms[i]
            typ = atm.atomtype
            if atm.reg_dim_atm[0] == []:
                if len(typ) == 1:
                    if typ[0].label == 'R':
                        extents.extend(self.specify_atom_extensions(i, basename, R))  # specify types of atoms
                    elif typ[0].label == 'R!H':
                        extents.extend(self.specify_atom_extensions(i, basename, RnH))
                    elif typ[0].label == 'Rx':
                        extents.extend(self.specify_atom_extensions(i, basename, r))
                else:
                    extents.extend(self.specify_atom_extensions(i, basename, typ))
            else:
                if len(typ) == 1:
                    if typ[0].label == 'R':
                        extents.extend(
                            self.specify_atom_extensions(i, basename, atm.reg_dim_atm[0]))  # specify types of atoms
                    elif typ[0].label == 'R!H':
                        extents.extend(self.specify_atom_extensions(i, basename, list(set(atm.reg_dim_atm[0]) & set(r))))
                else:
                    extents.extend(self.specify_atom_extensions(i, basename, list(set(typ) & set(atm.reg_dim_atm[0]))))
            if atm.reg_dim_u == []:
                if len(atm.radical_electrons) != 1:
                    if len(atm.radical_electrons) == 0:
                        extents.extend(self.specify_unpaired_extensions(i, basename, r_un))
                    else:
                        extents.extend(self.specify_unpaired_extensions(i, basename, atm.radical_electrons))
            else:
                if len(atm.radical_electrons) != 1 and len(atm.reg_dim_u[0]) != 1:
                    if len(atm.radical_electrons) == 0:
                        extents.extend(self.specify_unpaired_extensions(i, basename, atm.reg_dim_u[0]))
                    else:
                        extents.extend(self.specify_unpaired_extensions(i, basename, list(
                            set(atm.radical_electrons) & set(atm.reg_dim_u[0]))))
            if atm.reg_dim_r[0] == [] and 'inRing' not in atm.props:
                extents.extend(self.specify_ring_extensions(i, basename))

            extents.extend(self.specify_external_new_bond_extensions(i, basename, r_bonds))
            for j, atm2 in enumerate(atoms):
                if j < i and not self.has_bond(atm, atm2):
                    extents.extend(self.specify_internal_new_bond_extensions(i, j, n_splits, basename, r_bonds))
                elif j < i:
                    bd = self.get_bond(atm, atm2)
                    if len(bd.order) > 1 and bd.reg_dim == []:
                        extents.extend(self.specify_bond_extensions(i, j, basename, bd.order))
                    elif len(bd.order) > 1 and len(bd.reg_dim[0]) > 1 and len(bd.reg_dim[0]) > len(bd.reg_dim[1]):
                        extents.extend(self.specify_bond_extensions(i, j, basename, bd.reg_dim[0]))

        else:
            raise ValueError('atm_ind must be defined if atm_ind2 is defined')

        return extents

    def specify_atom_extensions(self, i, basename, r):
        """
        generates extensions for specification of the type of atom defined by a given atomtype
        or set of atomtypes
        """
        cython.declare(grps=list, labelList=list, Rset=set, item=AtomType, grp=Group, grpc=Group, k=AtomType, p=str)

        grps = []
        Rset = set(r)

        #consider node splitting        
        for item in r:
            grp = deepcopy(self)
            grpc = deepcopy(self)
            old_atom_type = grp.atoms[i].atomtype
            grp.atoms[i].atomtype = [item]
            grpc.atoms[i].atomtype = list(Rset - {item})

            if len(grpc.atoms[i].atomtype) == 0:
                grpc = None

            if len(old_atom_type) > 1:
                labelList = []
                old_atom_type_str = ''
                for k in old_atom_type:
                    labelList.append(k.label)
                for p in sorted(labelList):
                    old_atom_type_str += p
            elif len(old_atom_type) == 0:
                old_atom_type_str = ""
            else:
                old_atom_type_str = old_atom_type[0].label

            grps.append(
                (grp, grpc, basename + '_' + str(i + 1) + old_atom_type_str + '->' + item.label, 'atomExt', (i,)))

        #generate an extension without node splitting
        if len(self.atoms[i].atomtype)>len(Rset):
            print('generating a non-splitting extension')
            if all(r in self.atoms[i].atomtype for r in Rset): 
                #that means even if we update the atomtype of the atom to the Rset, it will still be a specification
                grp = deepcopy(self)
                grp.atoms[i].atomtype = list(Rset)
                
                #rename
                old_atom_type = grp.atoms[i].atomtype

                if len(old_atom_type) > 1:
                    labelList = []
                    old_atom_type_str = ''
                    for k in old_atom_type:
                        labelList.append(k.label)
                    for p in sorted(labelList):
                        old_atom_type_str += p
                elif len(old_atom_type) == 0:
                    old_atom_type_str = ""
                else:
                    old_atom_type_str = old_atom_type[0].label

                grps.append(
                (grp, None, basename + '_' + str(i + 1) + old_atom_type_str + '->' + ''.join(r.label for r in Rset), 'atomExt', (i,)))
       
        return grps

    def specify_ring_extensions(self, i, basename):
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

        atom_type = grp.atoms[i].atomtype

        if len(atom_type) > 1:
            atom_type_str = ''
            for k in atom_type:
                label_list.append(k.label)
            for p in sorted(label_list):
                atom_type_str += p
        elif len(atom_type) == 0:
            atom_type_str = ""
        else:
            atom_type_str = atom_type[0].label

        grps.append((grp, grpc, basename + '_' + str(i + 1) + atom_type_str + '-inRing', 'ringExt', (i,)))

        return grps

    def specify_unpaired_extensions(self, i, basename, r_un):
        """
        generates extensions for specification of the number of electrons on a given atom
        """

        grps = []
        label_list = []

        Rset = set(r_un)
        for item in r_un:
            grp = deepcopy(self)
            grpc = deepcopy(self)
            grp.atoms[i].radical_electrons = [item]
            grpc.atoms[i].radical_electrons = list(Rset - {item})

            if len(grpc.atoms[i].radical_electrons) == 0:
                grpc = None

            atom_type = grp.atoms[i].atomtype

            if len(atom_type) > 1:
                atom_type_str = ''
                for k in atom_type:
                    label_list.append(k.label)
                for p in sorted(label_list):
                    atom_type_str += p
            elif len(atom_type) == 0:
                atom_type_str = ""
            else:
                atom_type_str = atom_type[0].label

            grps.append((grp, grpc, basename + '_' + str(i + 1) + atom_type_str + '-u' + str(item), 'elExt', (i,)))

        return grps

    def specify_internal_new_bond_extensions(self, i, j, n_splits, basename, r_bonds):
        """
        generates extensions for creation of a bond (of undefined order)
        between two atoms indexed i,j that already exist in the group and are unbonded
        """
        cython.declare(newgrp=Group)

        label_list = []

        newgrp = deepcopy(self)
        newgrp.add_bond(GroupBond(newgrp.atoms[i], newgrp.atoms[j], r_bonds))

        atom_type_i = newgrp.atoms[i].atomtype
        atom_type_j = newgrp.atoms[j].atomtype

        if len(atom_type_i) > 1:
            atom_type_i_str = ''
            for k in atom_type_i:
                label_list.append(k.label)
            for k in sorted(label_list):
                atom_type_i_str += k
        elif len(atom_type_i) == 0:
            atom_type_i_str = ""
        else:
            atom_type_i_str = atom_type_i[0].label
        if len(atom_type_j) > 1:
            atom_type_j_str = ''
            for k in atom_type_j:
                label_list.append(k.label)
            for p in sorted(label_list):
                atom_type_j_str += p
        elif len(atom_type_j) == 0:
            atom_type_j_str = ""
        else:
            atom_type_j_str = atom_type_j[0].label

        if len(newgrp.split()) < n_splits:  # if this formed a bond between two seperate groups in the
            return []
        else:
            return [(newgrp, None,
                     basename + '_Int-' + str(i + 1) + atom_type_i_str + '-' + str(j + 1) + atom_type_j_str,
                     'intNewBondExt', (i, j))]

    def specify_external_new_bond_extensions(self, i, basename, r_bonds):
        """
        generates extensions for the creation of a bond (of undefined order) between
        an atom and a new atom that is not H
        """
        cython.declare(ga=GroupAtom, newgrp=Group, j=int)

        label_list = []

        ga = GroupAtom([ATOMTYPES['R!H']])
        newgrp = deepcopy(self)
        newgrp.add_atom(ga)
        j = newgrp.atoms.index(ga)
        newgrp.add_bond(GroupBond(newgrp.atoms[i], newgrp.atoms[j], r_bonds))
        atom_type = newgrp.atoms[i].atomtype
        if len(atom_type) > 1:
            atom_type_str = ''
            for k in atom_type:
                label_list.append(k.label)
            for p in sorted(label_list):
                atom_type_str += p
        elif len(atom_type) == 0:
            atom_type_str = ""
        else:
            atom_type_str = atom_type[0].label

        return [(newgrp, None, basename + '_Ext-' + str(i + 1) + atom_type_str + '-R', 'extNewBondExt',
                 (len(newgrp.atoms) - 1,))]

    def specify_bond_extensions(self, i, j, basename, r_bonds):
        """
        generates extensions for the specification of bond order for a given bond
        """
        cython.declare(grps=list, label_list=list, Rbset=set, bd=float, grp=Group, grpc=Group)
        grps = []
        label_list = []
        Rbset = set(r_bonds)
        bdict = {1: '-', 2: '=', 3: '#', 1.5: '-=', 4: '$', 0.05: '..', 0: '--'}

        for bd in r_bonds:
            grp = deepcopy(self)
            grpc = deepcopy(self)
            grp.atoms[i].bonds[grp.atoms[j]].order = [bd]
            grp.atoms[j].bonds[grp.atoms[i]].order = [bd]
            grpc.atoms[i].bonds[grpc.atoms[j]].order = list(Rbset - {bd})
            grpc.atoms[j].bonds[grpc.atoms[i]].order = list(Rbset - {bd})

            if len(list(Rbset - {bd})) == 0:
                grpc = None

            atom_type_i = grp.atoms[i].atomtype
            atom_type_j = grp.atoms[j].atomtype

            if len(atom_type_i) > 1:
                atom_type_i_str = ''
                for k in atom_type_i:
                    label_list.append(k.label)
                for p in sorted(label_list):
                    atom_type_i_str += p
            elif len(atom_type_i) == 0:
                atom_type_i_str = ""
            else:
                atom_type_i_str = atom_type_i[0].label
            if len(atom_type_j) > 1:
                atom_type_j_str = ''
                for k in atom_type_j:
                    label_list.append(k.label)
                for p in sorted(label_list):
                    atom_type_j_str += p
            elif len(atom_type_j) == 0:
                atom_type_j_str = ""
            else:
                atom_type_j_str = atom_type_j[0].label

            b = None
            for v in bdict.keys():
                if abs(v - bd) < 1e-4:
                    b = bdict[v]


            grps.append((grp, grpc,
                        basename + '_Sp-' + str(i + 1) + atom_type_i_str + b + str(j + 1) + atom_type_j_str,
                         'bondExt', (i, j)))

        return grps

    def clear_labeled_atoms(self):
        """
        Remove the labels from all atoms in the molecular group.
        """
        cython.declare(atom=GroupAtom)
        for atom in self.vertices:
            atom.label = ''

    def contains_labeled_atom(self, label):
        """
        Return ``True`` if the group contains an atom with the label
        `label` and ``False`` otherwise.
        """
        cython.declare(atom=GroupAtom)
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def get_labeled_atoms(self, label):
        """
        Return the atom in the group that is labeled with the given `label`.
        Raises :class:`ValueError` if no atom in the group has that label.
        """
        cython.declare(atom=GroupAtom, alist=list)
        alist = [atom for atom in self.vertices if atom.label == label]
        if alist == []:
            raise ValueError('No atom in the functional group \n{1}\n has the label '
                             '"{0}".'.format(label, self.to_adjacency_list()))
        return alist

    def get_all_labeled_atoms(self):
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
            for atomtype in atom.atomtype:
                if match is None:
                    # This is the first type in the list, so check all elements
                    for element in allElements:
                        if atomtype.is_specific_case_of(ATOMTYPES[element]):
                            match = element
                            break
                else:
                    # We've already matched one atomtype, now confirm that the rest are the same
                    if not atomtype.is_specific_case_of(ATOMTYPES[match]):
                        same = False
                        break
            # If match is None, then the group is not a specific case of any element
            if match is not None and same:
                if match in element_count:
                    element_count[match] += 1
                else:
                    element_count[match] = 1

        return element_count

    def from_adjacency_list(self, adjlist, check_consistency=True):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from rmgpy.molecule.adjlist import from_adjacency_list
        self.vertices, multiplicity, self.metal, self.facet = from_adjacency_list(adjlist, group=True, check_consistency=check_consistency)
        if multiplicity is not None:
            self.multiplicity = multiplicity
        self.update()
        return self

    def to_adjacency_list(self, label=''):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from rmgpy.molecule.adjlist import to_adjacency_list
        return to_adjacency_list(self.vertices, multiplicity=self.multiplicity, metal=self.metal, facet=self.facet, label=label, group=True)

    def update_fingerprint(self):
        """
        Update the molecular fingerprint used to accelerate the subgraph
        isomorphism checks.
        """
        cython.declare(atom=GroupAtom)

        self.elementCount = self.get_element_count()
        self.radicalCount = 0
        for atom in self.vertices:
            if len(atom.radical_electrons) >= 1:
                self.radicalCount += atom.radical_electrons[0]

    def is_isomorphic(self, other, initial_map=None, generate_initial_map=False, save_order=False, strict=True):
        """
        Returns ``True`` if two graphs are isomorphic and ``False``
        otherwise. The `initial_map` attribute can be used to specify a required
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
        return Graph.is_isomorphic(self, other, initial_map, generate_initial_map, save_order=save_order)

    def find_isomorphism(self, other, initial_map=None, save_order=False, strict=True):
        """
        Returns ``True`` if `other` is isomorphic and ``False``
        otherwise, and the matching mapping. The `initial_map` attribute can be
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
        return Graph.find_isomorphism(self, other, initial_map, save_order=save_order)

    def is_subgraph_isomorphic(self, other, initial_map=None, generate_initial_map=False, save_order=False):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` if self is more specific than other.
        The `initial_map` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short, m1=str, m2=str)
        cython.declare(a=GroupAtom, L=list)
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))

        group = other

        if generate_initial_map:
            keys = []
            atms = []
            initial_map = dict()
            for atom in self.atoms:
                if atom.label and atom.label != '':
                    L = [a for a in other.atoms if a.label == atom.label]
                    if L == []:
                        return False
                    elif len(L) == 1:
                        initial_map[atom] = L[0]
                    else:
                        keys.append(atom)
                        atms.append(L)
            if atms:
                for atmlist in itertools.product(*atms):
                    if len(set(atmlist)) != len(atmlist):
                        # skip entries that map multiple graph atoms to the same subgraph atom
                        continue
                    for i, key in enumerate(keys):
                        initial_map[key] = atmlist[i]
                    if (self.is_mapping_valid(other, initial_map, equivalent=False) and
                            Graph.is_subgraph_isomorphic(self, other, initial_map, save_order=save_order)):
                        return True
                else:
                    return False
            else:
                if not self.is_mapping_valid(other, initial_map, equivalent=False):
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
        if self.metal:
            for m1 in self.metal:
                if group.metal:
                    for m2 in group.metal:
                        if m1 == m2: break
                    else:
                        return False
        else:
            if group.metal: return False
        if self.facet:
            for m1 in self.facet:
                if group.facet:
                    for m2 in group.facet:
                        if m1 == m2: break
                    else:
                        return False
        else:
            if group.facet: return False
        # Do the isomorphism comparison
        return Graph.is_subgraph_isomorphic(self, other, initial_map, save_order=save_order)

    def find_subgraph_isomorphisms(self, other, initial_map=None, save_order=False):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. In other words, return ``True`` is self is more specific than other.
        Also returns the lists all of valid mappings. The
        `initial_map` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        cython.declare(group=Group)
        cython.declare(mult1=cython.short, mult2=cython.short, m1=str, m2=str)

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
        if self.metal:
            for m1 in self.metal:
                if group.metal:
                    for m2 in group.metal:
                        if m1 == m2: break
                    else:
                        return []
        else:
            if group.metal:
                return []
        if self.facet:
            for m1 in self.facet:
                if group.facet:
                    for m2 in group.facet:
                        if m1 == m2: break
                    else:
                        return []
        else:
            if group.facet:
                return []

        # Do the isomorphism comparison
        return Graph.find_subgraph_isomorphisms(self, other, initial_map, save_order=save_order)

    def is_identical(self, other, save_order=False):
        """
        Returns ``True`` if `other` is identical and ``False`` otherwise.
        The function `is_isomorphic` respects wildcards, while this function
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
        if not self.is_subgraph_isomorphic(other, None, save_order=save_order):
            return False
        elif not other.is_subgraph_isomorphic(self, None, save_order=save_order):
            return False
        else:
            return True

    def is_aromatic_ring(self):
        """
        This method returns a boolean telling if the group has a 5 or 6 cyclic with
        benzene bonds exclusively
        """

        ring_size = len(self.atoms)
        if ring_size not in [5, 6]:
            return False
        for ring_atom in self.atoms:
            for bonded_atom, bond in ring_atom.edges.items():
                if bonded_atom in self.atoms:
                    if not bond.is_benzene():
                        return False
        return True

    def has_wildcards(self):
        """
        This function is a Group level wildcards checker.

        Returns a 'True' if any of the atoms in this group has wildcards.
        """

        for atom1 in self.atoms:
            if atom1.has_wildcards():
                return True

        return False

    def standardize_atomtype(self):
        """
        This function changes the atomtypes in a group if the atom must
        be a specific atomtype based on its bonds and valency.

        Currently only standardizes oxygen, carbon and sulfur ATOMTYPES

        We also only check when there is exactly one atomtype,
        one bondType, one radical setting.
        For any group where there are wildcards or multiple attributes,
        we cannot apply this check.

        In the case where the atomtype is ambiguous based on bonds
        and valency, this function will not change the type.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """
        modified = False

        # If this atom or any of its ligands has wild cards, then don't try to standardize
        if self.has_wildcards(): return modified

        # list of :class:AtomType which are elements with more sub-divided atomtypes beneath them
        specifics = [elementLabel for elementLabel in allElements if elementLabel not in nonSpecifics]
        for atom in self.atoms:
            claimed_atom_type = atom.atomtype[0]
            new_atom_type = None
            element = None
            # Ignore elements that do not have more than one atomtype
            if claimed_atom_type.label in nonSpecifics: continue
            for elementLabel in specifics:
                if (claimed_atom_type.label == elementLabel or
                        ATOMTYPES[claimed_atom_type.label] in ATOMTYPES[elementLabel].specific):
                    element = ATOMTYPES[elementLabel]
                    break

            # claimed_atom_type is not in one of the specified elements
            if not element:
                continue
            # Don't standardize atomtypes for nitrogen for now
            # The work on the nitrogen atomtypes is still incomplete
            elif element is ATOMTYPES['N']:
                continue

            group_features = get_features(atom, atom.bonds)

            bond_order = atom.get_total_bond_order()
            filled_valency = atom.radical_electrons[0] + bond_order

            # For an atomtype to be known for certain, the valency must be filled
            # within 1 of the total valency available
            if filled_valency >= PeriodicSystem.valence_electrons[self.symbol][element] - 1:
                for specific_atom_type in element.specific:
                    atomtype_feature_list = specific_atom_type.get_features()
                    for mol_feature, atomtype_feature in zip(group_features, atomtype_feature_list):
                        if atomtype_feature == []:
                            continue
                        elif mol_feature not in atomtype_feature:
                            break
                    else:
                        if specific_atom_type is ATOMTYPES['Oa'] or specific_atom_type is ATOMTYPES['Sa']:
                            if atom.lone_pairs == 3 or atom.radical_electrons == 2:
                                new_atom_type = specific_atom_type
                                break
                        else:
                            new_atom_type = specific_atom_type
                            break

            # set the new atom type if the algorithm found one
            if new_atom_type and new_atom_type is not claimed_atom_type:
                atom.atomtype[0] = new_atom_type
                modified = True

        return modified

    def create_and_connect_atom(self, atomtypes, connecting_atom, bond_orders):
        """
        This method creates an non-radical, uncharged, :class:GroupAtom with specified list of atomtypes and
        connects it to one atom of the group, 'connecting_atom'. This is useful for making sample atoms.

        Args:
            atomtypes: list of atomtype labels (strs)
            connecting_atom: :class:GroupAtom that is connected to the new benzene atom
            bond_orders: list of bond Orders connecting new_atom and connecting_atom

        Returns: the newly created atom
        """
        atomtypes = [ATOMTYPES[label] for label in atomtypes]  # turn into :class: atomtype instead of labels

        new_atom = GroupAtom(atomtype=atomtypes, radical_electrons=[0], charge=[], label='', lone_pairs=None)
        new_bond = GroupBond(connecting_atom, new_atom, order=bond_orders)
        self.add_atom(new_atom)
        self.add_bond(new_bond)
        return new_atom

    def add_explicit_ligands(self):
        """
        This function O2d/S2d ligand to CO or CS atomtypes if they are not already there.

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        modified = False

        atoms_to_add_to = []

        for index, atom in enumerate(self.atoms):
            claimed_atom_type = atom.atomtype[0]
            # Do not perform is this atom has wildCards
            if atom.has_wildcards():
                continue
            elif claimed_atom_type is ATOMTYPES['CO'] or claimed_atom_type is ATOMTYPES['CS']:
                for bond12 in atom.bonds.values():
                    if bond12.is_double():
                        break
                else:
                    atoms_to_add_to.append(index)

        for atomIndex in atoms_to_add_to:
            modified = True
            atomtypes = None
            if self.atoms[atomIndex].atomtype[0] is ATOMTYPES['CO']:
                atomtypes = ['O2d']
            elif self.atoms[atomIndex].atomtype[0] is ATOMTYPES['CS']:
                atomtypes = ['S2d']
            self.create_and_connect_atom(atomtypes, self.atoms[atomIndex], [2])

        return modified

    def standardize_group(self):
        """
        This function modifies groups to make them have a standard AdjList form.

        Currently it makes atomtypes as specific as possible and makes CO/CS atomtypes
        have explicit O2d/S2d ligands. Other functions can be added as necessary

        Returns a 'True' if the group was modified otherwise returns 'False'
        """

        # If viable then we apply current conventions:
        check_list = []
        check_list.append(self.standardize_atomtype())
        check_list.append(self.add_explicit_ligands())
        return any(check_list)

    def add_implicit_atoms_from_atomtype(self):
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
            atomtype_feature_list = atom1.atomtype[0].get_features()
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
                if bond12.is_single():
                    single += 1
                elif bond12.is_double():
                    if atom2.is_oxygen():
                        o_double += 1
                    elif atom2.is_sulfur():
                        s_double += 1
                    else:
                        # r_double is for double bonds NOT to oxygen or Sulfur
                        r_double += 1
                elif bond12.is_triple():
                    triple += 1
                elif bond12.is_quadruple():
                    quadruple += 1
                elif bond12.is_benzene():
                    benzene += 1

            while o_double < o_double_required[0]:
                o_double += 1
                new_atom = GroupAtom(atomtype=[ATOMTYPES['O']], radical_electrons=[0], charge=[], label='',
                                     lone_pairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[2])
                implicit_atoms[new_atom] = new_bond
            while s_double < s_double_required[0]:
                s_double += 1
                new_atom = GroupAtom(atomtype=[ATOMTYPES['S']], radical_electrons=[0], charge=[], label='',
                                     lone_pairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[2])
                implicit_atoms[new_atom] = new_bond
            while r_double < r_double_required[0] or r_double + o_double + s_double < all_double_required[0]:
                r_double += 1
                new_atom = GroupAtom(atomtype=[ATOMTYPES['C']], radical_electrons=[0], charge=[], label='',
                                     lone_pairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[2])
                implicit_atoms[new_atom] = new_bond
            while triple < triple_required[0]:
                triple += 1
                new_atom = GroupAtom(atomtype=[ATOMTYPES['C']], radical_electrons=[0], charge=[], label='',
                                     lone_pairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[3])
                implicit_atoms[new_atom] = new_bond
            while quadruple < quadruple_required[0]:
                quadruple += 1
                new_atom = GroupAtom(atomtype=[ATOMTYPES['C']], radical_electrons=[0], charge=[], label='',
                                     lone_pairs=None)
                new_bond = GroupBond(atom1, new_atom, order=[4])
                implicit_atoms[new_atom] = new_bond

        for atom, bond in implicit_atoms.items():
            copy_group.add_atom(atom)
            copy_group.add_bond(bond)

        for atom, lone_pair in lone_pairs_required.items():
            if lone_pair: atom.lone_pairs = lone_pair

        return copy_group

    def classify_benzene_carbons(self, partners=None):
        """
        Args:
            group: :class:Group with atoms to classify
            partners: dictionary of partnered up atoms, which must be a cbf atom

        Some non-carbon 'benzene' type atoms (eg. N3b) are included and classified.

        Returns: tuple with lists of each atom classification:
        cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs
        """
        if not partners:
            partners = {}

        cb_atom_list = []
        cbf_atom_list = []  # All Cbf Atoms
        cbf_atom_list1 = []  # Cbf Atoms that are bonded to exactly one other Cbf (part of 2 rings)
        cbf_atom_list2 = []  # Cbf that are sandwiched between two other Cbf (part of 2 rings)
        connected_cbfs = {}  # dictionary of connections to other cbf Atoms

        # Only want to work with benzene bonds on carbon
        labels_of_carbon_atom_types = [x.label for x in ATOMTYPES['C'].specific] + ['C']
        # Also allow with R!H and some other aromatic groups
        labels_of_carbon_atom_types.extend(['R!H', 'N5b', 'N3b', 'N5bd', 'O4b', 'P3b', 'P5b', 'P5bd', 'S4b'])
        # Why are Sib and Sibf missing?

        for atom in self.atoms:
            atomtype = atom.atomtype[0]
            if atom.atomtype[0].label not in labels_of_carbon_atom_types:
                continue
            elif atom.atomtype[0].label in ['Cb', 'N5b', 'N3b', 'N5bd', 'O4b', 'P3b', 'P5b', 'P5bd', 'S4b']:  # Make Cb and N3b into normal cb atoms
                cb_atom_list.append(atom)
            elif atom.atomtype[0].label == 'Cbf':
                cbf_atom_list.append(atom)
            else:
                benzene_bonds = 0
                for atom2, bond12 in atom.bonds.items():
                    if bond12.is_benzene(): benzene_bonds += 1
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
                assert bond12.is_benzene(), "Cbf atom in {0} has a bond with an order other than 1.5".format(self)

        return cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs

    def add_implicit_benzene(self):
        """
        Returns: A modified group with any implicit benzene rings added

        This method currently does not if there are wildcards in atomtypes or bond orders
        The current algorithm also requires that all Cb and Cbf are atomtyped

        There are other cases where the algorithm doesn't work. For example whenever there
        are many dangling Cb or Cbf atoms not in a ring, it is likely fail. In the database test
        (the only use thus far), we will require that any group with more than 3 Cbfs have
        complete rings. This is much stricter than this method can handle, but right now
        this method cannot handle very general cases, so it is better to be conservative.

        Note that it also works on other aromatic atomtypes like N5bd etc.
        """
        # Note that atomtypes like N5bd are mostly referred to as Cb in this code,
        # which was first written for just carbon.

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
        (cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs) = copy_group.classify_benzene_carbons()

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
                    elif bond12.is_benzene():
                        bonds_are_benzene = [True if bond23.is_benzene() else False for bond23 in atom2.bonds.values()]
                        if all(bonds_are_benzene) and 0 in atom2.radical_electrons:
                            potential_partner = atom2
                # Make a Cb atom the partner, now marking it as a Cbfatom
                if potential_partner:
                    partners[cbf_atom] = potential_partner
                    partners[potential_partner] = cbf_atom
                # otherwise create a new atom to be the partner
                else:
                    new_atom = copy_group.create_and_connect_atom(['Cbf'], cbf_atom, [1.5])
                    partners[cbf_atom] = new_atom
                    partners[new_atom] = cbf_atom

        # reclassify all atoms since we may have added new ones
        cb_atom_list, cbf_atom_list, cbf_atom_list1, cbf_atom_list2, connected_cbfs = copy_group.classify_benzene_carbons(partners)

        """
        Step 3. Sort all lists by connectivity

        In the coming steps, we will sort Cb/Cbf atom into their benzene rings. If we cannot
        find a ring to sort an atom into, we will create a new ring containing that atom.
        It is important that we always check atoms that are already connected to existing rings
        before completely disconnected atoms. Otherwise, we will erroneously create new rings.
        """
        cb_atom_list = copy_group.sort_by_connectivity(cb_atom_list)
        cbf_atom_list1 = copy_group.sort_by_connectivity(cbf_atom_list1)
        cbf_atom_list2 = copy_group.sort_by_connectivity(cbf_atom_list2)

        """
        Step 4. Initalize the list of rings with any benzene rings that are already explicitly stated

        The variable rings is a list of lists. Each list in rings represents one full benzene rings,
        so it will eventually have six benzene carbons in it. Each ring's list will have the atoms
        sorted by connectivity, such that any atom is bonded to the atoms preceding and following it
        in the list. The first and last atom of the list will also be bonded together.

        """
        rings = [cycle for cycle in copy_group.get_all_cycles_of_size(6) if Group(atoms=cycle).is_aromatic_ring()]

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
                new_atom = copy_group.create_and_connect_atom(['Cb'], cbf_atom, [1.5])
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
                    if check_set(ring2, ring1): break
                    # Next try to merge the ringseed into rings
                    merge_ring = merge_overlapping_benzene_rings(ring2, ring1, 2)
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
                if check_set(ring, new_ring_seed):
                    in_ring += 1
            # move on to next cbf_atom if we found two rings
            if in_ring == 2: continue
            # try to merge into existing rings, if cbf1 is connected
            for index, ring in enumerate(rings):
                merge_ring_dict = {}
                merge_ring = merge_overlapping_benzene_rings(ring, new_ring_seed, 1)
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
                if check_set(ring, [cb_atom]):
                    in_ring += 1
            # move on to next ring cb_atom if in a ring
            if in_ring == 1: continue
            # check to see if can be merged to an existing ring
            for index, ring in enumerate(rings):
                merge_ring_dict = {}
                merge_ring = add_cb_atom_to_ring(ring, cb_atom)
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
                new_atom = copy_group.create_and_connect_atom(['Cb'], last_atom, [1.5])
                merged_ring_dict[index].append(new_atom)
                # At the end attach to the other endpoint
                if x == carbons_to_grow - 1:
                    new_bond = GroupBond(ring[0], new_atom, order=[1.5])
                    copy_group.add_bond(new_bond)

        return copy_group

    def pick_wildcards(self):
        """
        Returns: the :class:Group object without wildcards in either atomtype or bonding

        This function will naively pick the first atomtype for each atom, but will try
        to pick bond orders that make sense given the selected atomtypes
        """
        for atom1 in self.atoms:
            atom1.atomtype = [atom1.atomtype[0]]
            for atom2, bond12 in atom1.bonds.items():
                # skip dynamic bond ordering if there are no wildcards
                if len(bond12.order) < 2:
                    continue
                atom1_features = atom1.atomtype[0].get_features()
                atom2_features = atom2.atomtype[0].get_features()
                # skip dynamic bond ordering if there are no features required by the atomtype
                if not any(atom1_features) and not any(atom2_features):
                    bond12.order = [bond12.order[0]]
                    atom2.bonds[atom1].order = bond12.order
                    continue

                atom1_bonds = atom1.count_bonds()  # count_bonds list must match get_features list
                atom2_bonds = atom2.count_bonds()
                required_features1 = [atom1_features[x][0] - atom1_bonds[x] if atom1_features[x] else 0 for x in
                                      range(len(atom1_bonds))]
                required_features2 = [atom2_features[x][0] - atom2_bonds[x] if atom2_features[x] else 0 for x in
                                      range(len(atom2_bonds))]

                # subtract 1 from all_double for each s_double, o_double, r_double so that we don't count them twice
                required_features1[1] = required_features1[1] - required_features1[2] - required_features1[3] - \
                                        required_features1[4]
                required_features2[1] = required_features2[1] - required_features2[2] - required_features2[3] - \
                                        required_features2[4]
                # reverse it because coincidentally the reverse has good priority on what to check first
                required_features1.reverse()
                required_features2.reverse()

                # required features are a now list of [benzene, quadruple, triple, s_double, o_double, r_double, all_double, single]
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
                            if index == 3:  # s_double bonds
                                if (feature1 > 0 and atom2.is_sulfur()) or (feature2 > 0 and atom1.is_sulfur()):
                                    bond12.order = [2]
                                    atom2.bonds[atom1].order = bond12.order
                                    break
                            elif index == 4:  # oDoubleBonds
                                if (feature1 > 0 and atom2.is_oxygen()) or (feature2 > 0 and atom1.is_oxygen()):
                                    bond12.order = [2]
                                    atom2.bonds[atom1].order = bond12.order
                                    break
                            else:  # r_double or all_double necessary
                                bond12.order = [2]
                                atom2.bonds[atom1].order = bond12.order
                                break
                else:  # no features required, then pick the first order
                    bond12.order = [bond12.order[0]]
                    atom2.bonds[atom1].order = bond12.order

        # if we have wildcard atomtypes pick one based on ordering of allElements
        for atom in self.atoms:
            for elementLabel in allElements:
                if ATOMTYPES[elementLabel] in atom.atomtype[0].specific:
                    atom.atomtype = [ATOMTYPES[elementLabel]]
                    break

    def make_sample_molecule(self):
        """
        Returns: A sample class :Molecule: from the group
        """

        modified_group = self.copy(deep=True)

        # Remove all wildcards
        modified_group.pick_wildcards()

        # check that there are less than three Cbf atoms
        cbf_count = 0
        for atom in modified_group.atoms:
            if atom.atomtype[0] is ATOMTYPES['Cbf']: cbf_count += 1
        if cbf_count > 3:
            if not modified_group.is_benzene_explicit():
                raise ImplicitBenzeneError("{0} has more than three Cbf atoms and does not have fully explicit "
                                           "benzene rings.")

        # Add implicit atoms
        modified_group = modified_group.add_implicit_atoms_from_atomtype()

        # Add implicit benzene rings
        if not modified_group.is_benzene_explicit():
            modified_group = modified_group.add_implicit_benzene()
        # Make dictionary of :GroupAtoms: to :Atoms: and vice versa
        group_to_mol = {}
        mol_to_group = {}
        for atom in modified_group.atoms:
            mol_atom = atom.make_sample_atom()
            group_to_mol[atom] = mol_atom
            mol_to_group[mol_atom] = atom

        # create the molecule
        new_molecule = mol.Molecule(atoms=list(group_to_mol.values()))

        # Add explicit bonds to :Atoms:
        for atom1 in modified_group.atoms:
            for atom2, bond12 in atom1.bonds.items():
                bond12.make_bond(new_molecule, group_to_mol[atom1], group_to_mol[atom2])

        # Saturate up to expected valency
        for mol_atom in new_molecule.atoms:
            if mol_atom.charge:
                stated_charge = mol_atom.charge
            # otherwise assume no charge (or implicit atoms we assume hvae no charge)
            else:
                stated_charge = 0
            mol_atom.update_charge()
            if mol_atom.charge - stated_charge:
                hydrogen_needed = mol_atom.charge - stated_charge
                if mol_atom in mol_to_group and mol_to_group[mol_atom].atomtype[0].single:
                    max_single = max(mol_to_group[mol_atom].atomtype[0].single)
                    single_present = sum([1 for atom in mol_atom.bonds if mol_atom.bonds[atom].is_single()])
                    max_hydrogen = max_single - single_present
                    if hydrogen_needed > max_hydrogen: hydrogen_needed = max_hydrogen
                for x in range(hydrogen_needed):
                    new_h = mol.Atom('H', radical_electrons=0, lone_pairs=0, charge=0)
                    new_bond = mol.Bond(mol_atom, new_h, 1)
                    new_molecule.add_atom(new_h)
                    new_molecule.add_bond(new_bond)
                mol_atom.update_charge()

        new_molecule.update()

        # Check that the charge of atoms is expected
        for atom in new_molecule.atoms:
            if atom.charge != 0:
                if atom in mol_to_group:
                    group_atom = mol_to_group[atom]
                else:
                    raise UnexpectedChargeError(graph=new_molecule)
                # check hardcoded atomtypes
                positive_charged = ['H+',
                                    'Csc', 'Cdc',
                                    'N3sc', 'N5sc', 'N5dc', 'N5ddc', 'N5tc', 'N5b',
                                    'O2sc', 'O4sc', 'O4dc', 'O4tc',
                                    'P5sc', 'P5dc', 'P5ddc', 'P5tc', 'P5b',
                                    'S2sc', 'S4sc', 'S4dc', 'S4tdc', 'S6sc', 'S6dc', 'S6tdc']
                negative_charged = ['e',
                                    'C2sc', 'C2dc', 'C2tc',
                                    'N0sc', 'N1sc', 'N1dc', 'N5dddc',
                                    'O0sc',
                                    'P0sc', 'P1sc', 'P1dc', 'P5sc',
                                    'S0sc', 'S2sc', 'S2dc', 'S2tc', 'S4sc', 'S4dc', 'S4tdc', 'S6sc', 'S6dc', 'S6tdc']
                if atom.charge > 0 and any([group_atom.atomtype[0] is ATOMTYPES[x] or ATOMTYPES[x].is_specific_case_of(group_atom.atomtype[0]) for x in positive_charged]):
                    pass
                elif atom.charge < 0 and any([group_atom.atomtype[0] is ATOMTYPES[x] or ATOMTYPES[x].is_specific_case_of(group_atom.atomtype[0]) for x in negative_charged]):
                    pass
                elif atom.charge in group_atom.atomtype[0].charge:
                    # declared charge in original group is same as new charge
                    pass
                else:
                    raise UnexpectedChargeError(graph=new_molecule)

        return new_molecule

    def is_benzene_explicit(self):
        """

        Returns: 'True' if all Cb, Cbf (and other 'b' atoms)
        are in completely explicitly stated benzene rings.

        Otherwise return 'False'

        """
        # classify atoms
        cb_atom_list = []
        for atom in self.atoms:
            if sum(atom.atomtype[0].benzene): # atomtype has at least one benzene bond
                cb_atom_list.append(atom)
            else: # there may be some undeclared (eg. a generic C atomtype, with a benzene bond)
                for atom2, bond12 in atom.bonds.items():
                    if bond12.is_benzene():
                        cb_atom_list.append(atom)
                        break # can stop checking bonds for this atom

        # get all explicit benzene rings
        rings = [cycle for cycle in self.get_all_cycles_of_size(6) if Group(atoms=cycle).is_aromatic_ring()]

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

    def merge_groups(self, other, keep_identical_labels=False):
        """
        This function takes `other` :class:Group object and returns a merged :class:Group object based
        on overlapping labeled atoms between self and other

        Currently assumes `other` can be merged at the closest labelled atom
        if keep_identical_labels=True merge_groups will not try to merge atoms with the same labels
        """
        labeled1 = self.get_all_labeled_atoms()
        labeled2 = other.get_all_labeled_atoms()
        overlapping_labels = [x for x in labeled1 if x in labeled2]

        # dictionary of key = original atoms, value = copy atoms for deep copies of the two groups
        self_dict = self.copy_and_map()
        other_dict = other.copy_and_map()

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
        if not keep_identical_labels:
            bonds_to_remove = []
            for label in overlapping_labels:
                old_atom_b = self.get_labeled_atoms(label)[0]
                for old_atom_a, old_bond_ab in old_atom_b.bonds.items():
                    if old_atom_a.label not in overlapping_labels:  # this is bond we need to transfer over
                        # find and record bondAB from new backbone for later removal
                        new_atom_a = self_dict[old_atom_a]
                        new_atom_b = self_dict[old_atom_b]
                        new_atom_c = merged_group.get_labeled_atoms(old_atom_b.label)[0]
                        for atom, new_bond_ab in new_atom_a.bonds.items():
                            if atom is new_atom_b:
                                bonds_to_remove.append(new_bond_ab)
                                break
                        # add bond between atomA and AtomC
                        new_bond_ac = GroupBond(new_atom_a, new_atom_c, order=old_bond_ab.order)
                        merged_group.add_bond(new_bond_ac)
            # remove bonds from merged_group
            for bond in bonds_to_remove:
                merged_group.remove_bond(bond)

        return merged_group

    def reset_ring_membership(self):
        """
        Resets ring membership information in the GroupAtom.props attribute.
        """
        cython.declare(atom=GroupAtom)

        for atom in self.atoms:
            if 'inRing' in atom.props:
                del atom.props['inRing']
