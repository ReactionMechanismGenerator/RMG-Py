#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2020 Prof. William H. Green (whgreen@mit.edu),           #
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
This module contains functionality for reading from and writing to the
adjacency list format used by Reaction Mechanism Generator (RMG).
"""

import logging
import re
import warnings

from rmgpy.exceptions import InvalidAdjacencyListError
from rmgpy.molecule.atomtype import get_atomtype
from rmgpy.molecule.element import get_element, PeriodicSystem
from rmgpy.molecule.group import GroupAtom, GroupBond

import logging
from copy import deepcopy

import numpy as np

import rmgpy.molecule.element as elements
import rmgpy.molecule.group as gr
from rmgpy.molecule.atomtype import AtomType, ATOMTYPES, get_atomtype, AtomTypeError
from rmgpy.molecule.element import bdes
from rmgpy.molecule.graph import Vertex, Edge, Graph, get_vertex_connectivity_value


################################################################################

def kekulize():
    return None


def find_shortest_path():
    return None


bond_orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}

globals().update({
    'bond_orders': bond_orders,
})


class Atom(Vertex):
    """
    An atom. The attributes are:

    ==================== =================== ====================================
    Attribute            Type                Description
    ==================== =================== ====================================
    `atomtype`           :class:`AtomType`   The :ref:`atom type <atom-types>`
    `element`            :class:`Element`    The chemical element the atom represents
    `radical_electrons`  ``short``           The number of radical electrons
    `charge`             ``short``           The formal charge of the atom
    `label`              ``str``             A string label that can be used to tag individual atoms
    `coords`             ``numpy array``     The (x,y,z) coordinates in Angstrom
    `lone_pairs`         ``short``           The number of lone electron pairs
    `id`                 ``int``             Number assignment for atom tracking purposes
    `bonds`              ``dict``            Dictionary of bond objects with keys being neighboring atoms
    `props`              ``dict``            Dictionary for storing additional atom properties
    `mass`               ``int``             atomic mass of element (read only)
    `number`             ``int``             atomic number of element (read only)
    `symbol`             ``str``             atomic symbol of element (read only)
    ==================== =================== ====================================

    Additionally, the ``mass``, ``number``, and ``symbol`` attributes of the
    atom's element can be read (but not written) directly from the atom object,
    e.g. ``atom.symbol`` instead of ``atom.element.symbol``.
    """

    def __init__(self, element=None, radical_electrons=0, charge=0, label='', lone_pairs=-100, coords=np.array([], dtype=np.float128),
                 id=-1, props=None):
        Vertex.__init__(self)
        if isinstance(element, str):
            self.element = elements.__dict__[element]
        else:
            self.element = element
        self.radical_electrons = radical_electrons
        self.charge = charge
        self.label = label
        self.atomtype = None
        self.lone_pairs = lone_pairs
        self.coords = coords
        self.id = id
        self.props = props or {}

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '{0}{1}{2}'.format(
            str(self.element),
            '.' * self.radical_electrons,
            '+' * self.charge if self.charge > 0 else '-' * -self.charge,
        )

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "<Atom '{0}'>".format(str(self))

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
            'atomtype': self.atomtype.label if self.atomtype else None,
            'lone_pairs': self.lone_pairs,
        }
        if self.element.isotope == -1:
            element2pickle = self.element.symbol
        else:
            element2pickle = self.element
        return (Atom, (element2pickle, self.radical_electrons, self.charge, self.label), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an object.
        """
        self.edges = d['edges']
        self.connectivity1 = d['connectivity1']
        self.connectivity2 = d['connectivity2']
        self.connectivity3 = d['connectivity3']
        self.sorting_label = d['sorting_label']
        self.atomtype = ATOMTYPES[d['atomtype']] if d['atomtype'] else None
        self.lone_pairs = d['lone_pairs']

    def __hash__(self):
        """
        Define a custom hash method to allow Atom objects to be used in dictionaries and sets.
        """
        return hash(('Atom', self.symbol))

    def __eq__(self, other):
        """Method to test equality of two Atom objects."""
        return self is other

    def __lt__(self, other):
        """Define less than comparison. For comparing against other Atom objects (e.g. when sorting)."""
        if isinstance(other, Atom):
            return self.sorting_key < other.sorting_key
        else:
            raise NotImplementedError('Cannot perform less than comparison between Atom and '
                                      '{0}.'.format(type(other).__name__))

    def __gt__(self, other):
        """Define greater than comparison. For comparing against other Atom objects (e.g. when sorting)."""
        if isinstance(other, Atom):
            return self.sorting_key > other.sorting_key
        else:
            raise NotImplementedError('Cannot perform greater than comparison between Atom and '
                                      '{0}.'.format(type(other).__name__))

    @property
    def mass(self):
        return self.element.mass

    @property
    def number(self):
        return self.element.number

    @property
    def symbol(self):
        return self.element.symbol

    @property
    def bonds(self):
        return self.edges

    @property
    def sorting_key(self):
        """Returns a sorting key for comparing Atom objects. Read-only"""
        return self.number, -get_vertex_connectivity_value(self), self.radical_electrons, self.lone_pairs, self.charge

    def equivalent(self, other, strict=True):
        """
        Return ``True`` if `other` is indistinguishable from this atom, or
        ``False`` otherwise. If `other` is an :class:`Atom` object, then all
        attributes except `label` and 'ID' must match exactly. If `other` is an
        :class:`GroupAtom` object, then the atom must match any of the
        combinations in the atom pattern. If ``strict`` is ``False``, then only
        the element is compared and electrons are ignored.
        """
        if isinstance(other, Atom):
            atom = other
            if strict:
                return (self.element is atom.element
                        and self.radical_electrons == atom.radical_electrons
                        and self.lone_pairs == atom.lone_pairs
                        and self.charge == atom.charge
                        and self.atomtype is atom.atomtype)
            else:
                return self.element is atom.element
        elif isinstance(other, gr.GroupAtom):
            if not strict:
                raise NotImplementedError('There is currently no implementation of '
                                          'the strict argument for Group objects.')
            ap = other
            for a in ap.atomtype:
                if self.atomtype.equivalent(a): break
            else:
                return False
            if ap.radical_electrons:
                for radical in ap.radical_electrons:
                    if self.radical_electrons == radical: break
                else:
                    return False
            if ap.lone_pairs:
                for lp in ap.lone_pairs:
                    if self.lone_pairs == lp: break
                else:
                    return False
            if ap.charge:
                for charge in ap.charge:
                    if self.charge == charge: break
                else:
                    return False
            if 'inRing' in self.props and 'inRing' in ap.props:
                if self.props['inRing'] != ap.props['inRing']:
                    return False
            return True

    def is_specific_case_of(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. If `other` is an :class:`Atom` object, then this is the same
        as the :meth:`equivalent()` method. If `other` is an
        :class:`GroupAtom` object, then the atom must match or be more
        specific than any of the combinations in the atom pattern.
        """
        if isinstance(other, Atom):
            return self.equivalent(other)
        elif isinstance(other, gr.GroupAtom):
            atom = other
            if self.atomtype is None:
                return False
            for a in atom.atomtype:
                if self.atomtype.is_specific_case_of(a):
                    break
            else:
                return False
            if atom.radical_electrons:
                for radical in atom.radical_electrons:
                    if self.radical_electrons == radical:
                        break
                else:
                    return False
            if atom.lone_pairs:
                for lp in atom.lone_pairs:
                    if self.lone_pairs == lp:
                        break
                else:
                    return False
            if atom.charge:
                for charge in atom.charge:
                    if self.charge == charge:
                        break
                else:
                    return False
            if 'inRing' in self.props and 'inRing' in atom.props:
                if self.props['inRing'] != atom.props['inRing']:
                    return False
            elif 'inRing' not in self.props and 'inRing' in atom.props:
                return False
            return True

    def copy(self):
        """
        Generate a deep copy of the current atom. Modifying the
        attributes of the copy will not affect the original.
        """
        # a = Atom(self.element, self.radical_electrons, self.spin_multiplicity, self.charge, self.label)
        a = Atom.__new__(Atom)
        a.edges = {}
        a.reset_connectivity_values()
        a.element = self.element
        a.radical_electrons = self.radical_electrons
        a.charge = self.charge
        a.label = self.label
        a.atomtype = self.atomtype
        a.lone_pairs = self.lone_pairs
        a.coords = self.coords[:]
        a.id = self.id
        a.props = deepcopy(self.props)
        return a

    def is_hydrogen(self):
        """
        Return ``True`` if the atom represents a hydrogen atom or ``False`` if
        not.
        """
        return self.element.number == 1

    def is_non_hydrogen(self):
        """
        Return ``True`` if the atom does not represent a hydrogen atom or
        ``False`` if it does.
        """
        return self.element.number != 1

    def is_halogen(self):
        """
        Return ``True`` if the atom represents a halogen atom (F, Cl, Br, I)
        ``False`` if it does.
        """
        return self.element.number in [9, 17, 35, 53]

    def is_carbon(self):
        """
        Return ``True`` if the atom represents a carbon atom or ``False`` if
        not.
        """
        return self.element.number == 6

    def is_nitrogen(self):
        """
        Return ``True`` if the atom represents a nitrogen atom or ``False`` if
        not.
        """
        return self.element.number == 7

    def is_oxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if
        not.
        """
        return self.element.number == 8

    def is_fluorine(self):
        """
        Return ``True`` if the atom represents a fluorine atom or ``False`` if
        not.
        """
        return self.element.number == 9

    def is_surface_site(self):
        """
        Return ``True`` if the atom represents a surface site or ``False`` if not.
        """
        return self.symbol == 'X'

    def is_silicon(self):
        """
        Return ``True`` if the atom represents a silicon atom or ``False`` if
        not.
        """
        return self.element.number == 14

    def is_phosphorus(self):
        """
        Return ``True`` if the atom represents a phosphorus atom or ``False`` if
        not.
        """
        return self.element.number == 15

    def is_sulfur(self):
        """
        Return ``True`` if the atom represents a sulfur atom or ``False`` if
        not.
        """
        return self.element.number == 16

    def is_chlorine(self):
        """
        Return ``True`` if the atom represents a chlorine atom or ``False`` if
        not.
        """
        return self.element.number == 17

    def is_bromine(self):
        """
        Return ``True`` if the atom represents a bromine atom or ``False`` if
        not.
        """
        return self.element.number == 35

    def is_iodine(self):
        """
        Return ``True`` if the atom represents an iodine atom or ``False`` if
        not.
        """
        return self.element.number == 53

    def is_nos(self):
        """
        Return ``True`` if the atom represent either nitrogen, sulfur, or oxygen
        ``False`` if it does not.
        """
        return self.element.number in [7, 8, 16]

    def increment_radical(self):
        """
        Update the atom pattern as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        # Set the new radical electron count
        self.radical_electrons += 1
        if self.radical_electrons <= 0:
            raise gr.ActionError('Unable to update Atom due to GAIN_RADICAL action: '
                                 'Invalid radical electron set "{0}".'.format(self.radical_electrons))

    def decrement_radical(self):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        # Set the new radical electron count
        radical_electrons = self.radical_electrons = self.radical_electrons - 1
        if radical_electrons < 0:
            raise gr.ActionError('Unable to update Atom due to LOSE_RADICAL action: '
                                 'Invalid radical electron set "{0}".'.format(self.radical_electrons))

    def set_lone_pairs(self, lone_pairs):
        """
        Set the number of lone electron pairs.
        """
        # Set the number of electron pairs
        self.lone_pairs = lone_pairs
        if self.lone_pairs < 0:
            raise gr.ActionError('Unable to update Atom due to set_lone_pairs: '
                                 'Invalid lone electron pairs set "{0}".'.format(self.set_lone_pairs))
        self.update_charge()

    def increment_lone_pairs(self):
        """
        Update the lone electron pairs pattern as a result of applying a GAIN_PAIR action.
        """
        # Set the new lone electron pairs count
        self.lone_pairs += 1
        if self.lone_pairs <= 0:
            raise gr.ActionError('Unable to update Atom due to GAIN_PAIR action: '
                                 'Invalid lone electron pairs set "{0}".'.format(self.lone_pairs))
        self.update_charge()

    def decrement_lone_pairs(self):
        """
        Update the lone electron pairs pattern as a result of applying a LOSE_PAIR action.
        """
        # Set the new lone electron pairs count
        self.lone_pairs -= 1
        if self.lone_pairs < 0:
            raise gr.ActionError('Unable to update Atom due to LOSE_PAIR action: '
                                 'Invalid lone electron pairs set "{0}".'.format(self.lone_pairs))
        self.update_charge()

    def update_charge(self):
        """
        Update self.charge, according to the valence, and the
        number and types of bonds, radicals, and lone pairs.
        """
        if self.is_surface_site():
            self.charge = 0
            return
        valence_electron = elements.PeriodicSystem.valence_electrons[self.symbol]
        order = self.get_total_bond_order()
        self.charge = valence_electron - order - self.radical_electrons - 2 * self.lone_pairs

    def apply_action(self, action):
        """
        Update the atom pattern as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        # Invalidate current atom type
        self.atomtype = None
        act = action[0].upper()
        # Modify attributes if necessary
        if act in ['CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND']:
            # Nothing else to do here
            pass
        elif act == 'GAIN_RADICAL':
            for i in range(action[2]): self.increment_radical()
        elif act == 'LOSE_RADICAL':
            for i in range(abs(action[2])): self.decrement_radical()
        elif action[0].upper() == 'GAIN_PAIR':
            for i in range(action[2]): self.increment_lone_pairs()
        elif action[0].upper() == 'LOSE_PAIR':
            for i in range(abs(action[2])): self.decrement_lone_pairs()
        else:
            raise gr.ActionError('Unable to update Atom: Invalid action {0}".'.format(action))

    def get_total_bond_order(self):
        """
        This helper function is to help calculate total bond orders for an
        input atom.

        Some special consideration for the order `B` bond. For atoms having
        three `B` bonds, the order for each is 4/3.0, while for atoms having other
        than three `B` bonds, the order for  each is 3/2.0
        """
        num_b_bond = 0
        order = 0
        for bond in self.bonds.values():
            if bond.is_benzene():
                num_b_bond += 1
            else:
                order += bond.order

        if num_b_bond == 3:
            order += num_b_bond * 4 / 3.0
        else:
            order += num_b_bond * 3 / 2.0

        return order


################################################################################

class Bond(Edge):
    """
    A chemical bond. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``float``             The :ref:`bond type <bond-types>`
    `atom1`             ``Atom``              An Atom object connecting to the bond
    `atom2`             ``Atom``              An Atom object connecting to the bond
    =================== =================== ====================================

    """

    def __init__(self, atom1, atom2, order=1):
        Edge.__init__(self, atom1, atom2)
        if isinstance(order, str):
            self.set_order_str(order)
        else:
            self.order = order

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return self.get_order_str()

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return '<Bond "{0}">'.format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Bond, (self.vertex1, self.vertex2, self.order))

    def __hash__(self):
        """
        Define a custom hash method to allow Bond objects to be used in dictionaries and sets.
        """
        return hash(('Bond', self.order,
                     self.atom1.symbol if self.atom1 is not None else '',
                     self.atom2.symbol if self.atom2 is not None else ''))

    def __eq__(self, other):
        """Method to test equality of two Bond objects."""
        return self is other

    def __lt__(self, other):
        """Define less than comparison. For comparing against other Bond objects (e.g. when sorting)."""
        if isinstance(other, Bond):
            return self.sorting_key < other.sorting_key
        else:
            raise NotImplementedError('Cannot perform less than comparison between Bond and '
                                      '{0}.'.format(type(other).__name__))

    def __gt__(self, other):
        """Define greater than comparison. For comparing against other Bond objects (e.g. when sorting)."""
        if isinstance(other, Bond):
            return self.sorting_key > other.sorting_key
        else:
            raise NotImplementedError('Cannot perform greater than comparison between Bond and '
                                      '{0}.'.format(type(other).__name__))

    @property
    def atom1(self):
        return self.vertex1

    @property
    def atom2(self):
        return self.vertex2

    @property
    def sorting_key(self):
        """Returns a sorting key for comparing Bond objects. Read-only"""
        return (self.order,
                self.atom1.number if self.atom1 is not None else 0,
                self.atom2.number if self.atom2 is not None else 0)

    def get_bde(self):
        """
        estimate the bond dissociation energy in J/mol of the bond based on the order of the bond
        and the atoms involved in the bond
        """
        try:
            return bdes[(self.atom1.element.symbol, self.atom2.element.symbol, self.order)]
        except KeyError:
            raise KeyError('Bond Dissociation energy not known for combination: '
                           '({0},{1},{2})'.format(self.atom1.element.symbol, self.atom2.element.symbol, self.order))

    def equivalent(self, other):
        """
        Return ``True`` if `other` is indistinguishable from this bond, or
        ``False`` otherwise. `other` can be either a :class:`Bond` or a
        :class:`GroupBond` object.
        """
        if isinstance(other, Bond):
            bond = other
            return self.is_order(bond.get_order_num())
        elif isinstance(other, gr.GroupBond):
            bp = other
            return any([self.is_order(otherOrder) for otherOrder in bp.get_order_num()])

    def is_specific_case_of(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. `other` can be either a :class:`Bond` or a
        :class:`GroupBond` object.
        """
        # There are no generic bond types, so is_specific_case_of is the same as equivalent
        return self.equivalent(other)

    def get_order_str(self):
        """
        returns a string representing the bond order
        """
        if self.is_single():
            return 'S'
        elif self.is_benzene():
            return 'B'
        elif self.is_double():
            return 'D'
        elif self.is_triple():
            return 'T'
        elif self.is_quadruple():
            return 'Q'
        elif self.is_van_der_waals():
            return 'vdW'
        elif self.is_hydrogen_bond():
            return 'H'
        else:
            raise ValueError("Bond order {} does not have string representation.".format(self.order))

    def set_order_str(self, new_order):
        """
        set the bond order using a valid bond-order character
        """
        if new_order == 'S':
            self.order = 1
        elif new_order == 'D':
            self.order = 2
        elif new_order == 'T':
            self.order = 3
        elif new_order == 'B':
            self.order = 1.5
        elif new_order == 'Q':
            self.order = 4
        elif new_order == 'vdW':
            self.order = 0
        elif new_order == 'H':
            self.order = 0.1
        else:
            # try to see if an float disguised as a string was input by mistake
            try:
                self.order = float(new_order)
            except ValueError:
                raise TypeError('Bond order {} is not hardcoded into this method'.format(new_order))

    def get_order_num(self):
        """
        returns the bond order as a number
        """

        return self.order

    def set_order_num(self, new_order):
        """
        change the bond order with a number
        """

        self.order = new_order

    def copy(self):
        """
        Generate a deep copy of the current bond. Modifying the
        attributes of the copy will not affect the original.
        """
        # return Bond(self.vertex1, self.vertex2, self.order)
        b = Bond.__new__(Bond)
        b.vertex1 = self.vertex1
        b.vertex2 = self.vertex2
        b.order = self.order
        return b

    def is_van_der_waals(self):
        """
        Return ``True`` if the bond represents a van der Waals bond or
        ``False`` if not.
        """
        return self.is_order(0)

    def is_order(self, other_order):
        """
        Return ``True`` if the bond is of order other_order or ``False`` if
        not. This compares floats that takes into account floating point error

        NOTE: we can replace the absolute value relation with math.isclose when
        we swtich to python 3.5+
        """
        return abs(self.order - other_order) <= 1e-4

    def is_single(self):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not.
        """
        return self.is_order(1)

    def is_double(self):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not.
        """
        return self.is_order(2)

    def is_triple(self):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not.
        """
        return self.is_order(3)

    def is_quadruple(self):
        """
        Return ``True`` if the bond represents a quadruple bond or ``False`` if
        not.
        """
        return self.is_order(4)

    def is_benzene(self):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not.
        """
        return self.is_order(1.5)

    def is_hydrogen_bond(self):
        """
        Return ``True`` if the bond represents a hydrogen bond or ``False`` if
        not.
        """
        return self.is_order(0.1)

    def increment_order(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        increase the order by one.
        """
        if self.order <= 3.0001:
            self.order += 1
        else:
            raise gr.ActionError('Unable to increment Bond due to CHANGE_BOND action: '
                                 'Bond order "{0}" is greater than 3.'.format(self.order))

    def decrement_order(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        decrease the order by one.
        """
        if self.order >= 0.9999:
            self.order -= 1
        else:
            raise gr.ActionError('Unable to decrease Bond due to CHANGE_BOND action: '
                                 'bond order "{0}" is less than 1.'.format(self.order))

    def _change_bond(self, order):
        """
        Update the bond as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and can be any real number.
        """
        self.order += order
        if self.order < -0.0001 or self.order > 4.0001:
            raise gr.ActionError('Unable to update Bond due to CHANGE_BOND action: '
                                 'Invalid resulting order "{0}".'.format(self.order))

    def apply_action(self, action):
        """
        Update the bond as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            if isinstance(action[2], str):
                self.set_order_str(action[2])
            else:
                try:  # try to see if addable
                    self._change_bond(action[2])
                except TypeError:
                    raise gr.ActionError('Unable to update Bond due to CHANGE_BOND action: '
                                         'Invalid order "{0}".'.format(action[2]))
        else:
            raise gr.ActionError('Unable to update GroupBond: Invalid action {0}.'.format(action))

    def get_bond_string(self):
        """
        Represent the bond object as a string (eg. 'C#N'). The returned string is independent of the atom ordering, with
        the atom labels in alphabetical order (i.e. 'C-H' is possible but not 'H-C')
        :return: str
        """
        bond_symbol_mapping = {0.1: '~', 1: '-', 1.5: ':', 2: '=', 3: '#'}
        atom_labels = [self.atom1.symbol, self.atom2.symbol]
        atom_labels.sort()
        try:
            bond_symbol = bond_symbol_mapping[self.get_order_num()]
        except KeyError:
            # Direct lookup didn't work, but before giving up try
            # with the is_order() method which allows a little latitude
            # for floating point errors.
            for order, symbol in bond_symbol_mapping.items():
                if self.is_order(order):
                    bond_symbol = symbol
                    break
            else:  # didn't break
                bond_symbol = '<bond order {0}>'.format(self.get_order_num())
        return '{0}{1}{2}'.format(atom_labels[0], bond_symbol, atom_labels[1])


#################################################################################


class Saturator(object):
    @staticmethod
    def saturate(atoms):
        """
        Returns a list of atoms that is extended
        (and bond attributes) by saturating the valency of the non-hydrogen atoms with an
        appropriate number of hydrogen atoms.

        The required number of hydrogen atoms per heavy atom is determined as follows:
        H's =     max number of valence electrons - atom.radical_electrons
                    - 2* atom.lone_pairs - order - atom.charge

        """
        new_atoms = []
        for atom in atoms:
            try:
                max_number_of_valence_electrons = PeriodicSystem.valence_electrons[atom.symbol]
            except KeyError:
                raise InvalidAdjacencyListError(
                    'Cannot add hydrogens to adjacency list: Unknown orbital for atom "{0}".'.format(atom.symbol))

            order = atom.get_total_bond_order()

            number_of_h_to_be_added = max_number_of_valence_electrons - atom.radical_electrons - 2 * atom.lone_pairs - int(
                order) - atom.charge

            if number_of_h_to_be_added < 0:
                raise InvalidAdjacencyListError('Incorrect electron configuration on atom.')

            for _ in range(number_of_h_to_be_added):
                a = Atom(element='H', radical_electrons=0, charge=0, label='', lone_pairs=0)
                b = Bond(atom, a, 'S')
                new_atoms.append(a)
                atom.bonds[a] = b
                a.bonds[atom] = b
        atoms.extend(new_atoms)


class ConsistencyChecker(object):

    @staticmethod
    def check_partial_charge(atom):
        """
        Checks whether the partial charge attribute of the atom checks out with
        the theoretical one:

        """
        if atom.symbol == 'X':
            return  # because we can't check it.

        valence = PeriodicSystem.valence_electrons[atom.symbol]
        order = atom.get_total_bond_order()

        theoretical = valence - order - atom.radical_electrons - 2 * atom.lone_pairs

        if not (-0.301 < atom.charge - theoretical < 0.301):
            # It should be 0, but -0.1 is caused by a Hydrogen bond
            raise InvalidAdjacencyListError(
                'Invalid valency for atom {symbol} ({type}) with {radicals} unpaired electrons, '
                '{lone_pairs} pairs of electrons, {charge} charge, and bonds [{bonds}].'.format(
                    symbol=atom.symbol,
                    type=get_atomtype(atom, atom.edges).label,
                    radicals=atom.radical_electrons,
                    lone_pairs=atom.lone_pairs,
                    charge=atom.charge,
                    bonds=','.join([str(bond.order) for bond in atom.bonds.values()])
                )
            )

    @staticmethod
    def check_multiplicity(n_rad, multiplicity):
        """
        Check that the multiplicity complies with the formula: m = 2s + 1,
        where s is the sum of the spin [+/- (1/2) ] of the unpaired electrons

        For a simple radical (n_rad = 1):
        s = +1/2 , m = 2 (doublet)

        For a biradical, s can be either 0 [+0.5 + (-0.5) ] or 1 [+0.5 + (+0.5) ]
        and m = 1 (singlet) or m = 3 (triplet).
        """
        if n_rad in [0, 1]:
            if multiplicity != (n_rad + 1):
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        elif n_rad == 2:
            if not int(multiplicity) in [1, 3]:
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        elif n_rad == 3:
            if not int(multiplicity) in [4, 2]:
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        elif n_rad == 4:
            if not int(multiplicity) in [5, 3, 1]:
                raise InvalidAdjacencyListError('Multiplicity {0} not in agreement with total number of '
                                                'radicals {1}.'.format(multiplicity, n_rad))
        else:
            logging.warning("Consistency checking of multiplicity of molecules with "
                            "more than 4 unpaired electrons is not implemented yet!")

    @staticmethod
    def check_hund_rule(atom, multiplicity):
        """
        It is checked whether atoms with 2 unpaired electrons on the same atom
        result in a multiplicity of 3, and not 1.

        Unpaired electrons in 2 different orbitals belonging to the same atom
        should have the same spin, and hence, should result in a multiplicity of 3.
        """
        if atom.radical_electrons == 2 and multiplicity == 1:
            raise InvalidAdjacencyListError(
                "Violation of hund's rule. Invalid multiplicity of {0} because there is an "
                "atom with {1} unpaired electrons".format(multiplicity, atom.radical_electrons))


################################################################################

def from_old_adjacency_list(adjlist, group=False, saturate_h=False):
    """
    Convert a pre-June-2014 string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects. 
    It can read both "old style" that existed for years, an the "intermediate style" that
    existed for a few months in 2014, with the extra column of integers for lone pairs.
    """
    atoms = []
    atomdict = {}
    bonds = {}

    try:
        adjlist = adjlist.strip()
        lines = adjlist.splitlines()
        if adjlist == '' or len(lines) == 0:
            raise InvalidAdjacencyListError('Empty adjacency list.')

        # Skip the first line if it contains a label
        if len(lines[0].split()) == 1:
            label = lines.pop(0)
            if len(lines) == 0:
                raise InvalidAdjacencyListError("""Error in adjacency list\n{0}\nNo atoms specified.""".format(adjlist))

        mistake1 = re.compile(r'\{[^}]*\s+[^}]*\}')
        atomic_multiplicities = {}  # these are no longer stored on atoms, so we make a separate dictionary
        # Iterate over the remaining lines, generating Atom or GroupAtom objects
        for line in lines:

            # Sometimes people put spaces after commas, which messes up the
            # parse-by-whitespace. Examples include '{Cd, Ct}'.
            if mistake1.search(line):
                raise InvalidAdjacencyListError("Error in adjacency list: \n{1}\nspecies shouldn't have spaces inside "
                                                "braces: {0}".format(mistake1.search(line).group(), adjlist))

            # Sometimes commas are used to delimit bonds in the bond list,
            # so replace them just in case
            line = line.replace('},{', '} {')

            data = line.split()

            # Skip if blank line
            if len(data) == 0:
                continue

            # First item is index for atom
            # Sometimes these have a trailing period (as if in a numbered list),
            # so remove it just in case
            aid = int(data[0].strip('.'))

            # If second item starts with '*', then atom is labeled
            label = ''
            index = 1
            if data[1][0] == '*':
                label = data[1]
                index += 1

            # Next is the element or functional group element
            # A list can be specified with the {,} syntax
            atom_type = data[index]
            if atom_type[0] == '{':
                atom_type = atom_type[1:-1].split(',')
            else:
                atom_type = [atom_type]
            index += 1

            # Next is the electron state
            radical_electrons = []
            additional_lone_pairs = []
            elec_state = data[index].upper()
            if elec_state[0] == '{':
                elec_state = elec_state[1:-1].split(',')
            else:
                elec_state = [elec_state]
            if len(elec_state) == 0:
                raise InvalidAdjacencyListError(
                    "Error in adjacency list:\n{0}\nThere must be some electronic state defined for an "
                    "old adjlist".format(adjlist))
            for e in elec_state:
                if e == '0':
                    radical_electrons.append(0)
                    additional_lone_pairs.append(0)
                elif e == '1':
                    radical_electrons.append(1)
                    additional_lone_pairs.append(0)
                elif e == '2': 
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = 2 is not specific enough. "
                            "Please use 2S or 2T.".format(adjlist))
                    # includes 2S and 2T
                    radical_electrons.append(0); additional_lone_pairs.append(1)
                    radical_electrons.append(2); additional_lone_pairs.append(0)
                elif e == '2S':
                    radical_electrons.append(0); additional_lone_pairs.append(1)
                elif e == '2T':
                    radical_electrons.append(2); additional_lone_pairs.append(0)
                elif e == '3':
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = 3 is not specific enough. "
                            "Please use 3D or 3Q.".format(adjlist))
                    # includes 3D and 3Q
                    radical_electrons.append(1); additional_lone_pairs.append(1)
                    radical_electrons.append(3); additional_lone_pairs.append(0)
                elif e == '3D':
                    radical_electrons.append(1); additional_lone_pairs.append(1)
                elif e == '3Q':
                    radical_electrons.append(3); additional_lone_pairs.append(0)
                elif e == '4':
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = 4 is not specific enough. "
                            "Please use 4S, 4T, or 4V.".format(adjlist))
                    # includes 4S, 4T, and 4V
                    radical_electrons.append(0); additional_lone_pairs.append(2)
                    radical_electrons.append(2); additional_lone_pairs.append(1)
                    radical_electrons.append(4); additional_lone_pairs.append(0)
                elif e == '4S':
                    radical_electrons.append(0); additional_lone_pairs.append(2)
                elif e == '4T':
                    radical_electrons.append(2); additional_lone_pairs.append(1)
                elif e == '4V':
                    radical_electrons.append(4); additional_lone_pairs.append(0)
                elif e == 'X':
                    if not group:
                        raise InvalidAdjacencyListError(
                            "Error in adjacency list:\n{0}\nNumber of radical electrons = X is not specific enough. "
                            "Wildcards should only be used for groups.".format(adjlist))
                    radical_electrons = []
            index += 1

            # Next number defines the number of lone electron pairs (if provided)
            lone_pairs_of_electrons = None
            if len(data) > index:
                lp_state = data[index]
                if lp_state[0] == '{':
                    # this is the start of the chemical bonds - no lone pair info was provided
                    lone_pairs_of_electrons = None
                else:
                    if lp_state == '0':
                        lone_pairs_of_electrons = 0
                    if lp_state == '1':
                        lone_pairs_of_electrons = 1
                    if lp_state == '2':
                        lone_pairs_of_electrons = 2
                    if lp_state == '3':
                        lone_pairs_of_electrons = 3
                    if lp_state == '4':
                        lone_pairs_of_electrons = 4
                    index += 1
            else:  # no bonds or lone pair info provided.
                lone_pairs_of_electrons = None

            # Create a new atom based on the above information
            if group:
                if lone_pairs_of_electrons is not None:
                    lone_pairs_of_electrons = [additional + lone_pairs_of_electrons for additional in additional_lone_pairs]
                atom = GroupAtom(atomtype=atom_type,
                                 radical_electrons=sorted(set(radical_electrons)),
                                 charge=None,
                                 label=label,
                                 lone_pairs=lone_pairs_of_electrons,
                                 # Assign lone_pairs_of_electrons as None if it is not explicitly provided
                                 )

            else:
                if lone_pairs_of_electrons is not None:
                    # Intermediate adjlist representation
                    lone_pairs_of_electrons = lone_pairs_of_electrons + additional_lone_pairs[0]
                else:
                    # Add the standard number of lone pairs with the additional lone pairs
                    lone_pairs_of_electrons = PeriodicSystem.lone_pairs[atom_type[0]] + additional_lone_pairs[0]

                atom = Atom(element=atom_type[0],
                            radical_electrons=radical_electrons[0],
                            charge=0,
                            label=label,
                            lone_pairs=lone_pairs_of_electrons,
                            )
            # Add the atom to the list
            atoms.append(atom)
            atomdict[aid] = atom

            # Process list of bonds
            bonds[aid] = {}
            for datum in data[index:]:

                # Sometimes commas are used to delimit bonds in the bond list,
                # so strip them just in case
                datum = datum.strip(',')

                aid2, comma, order = datum[1:-1].partition(',')
                aid2 = int(aid2)
                if aid == aid2:
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAttempted to create a bond between '
                                                    'atom {0:d} and itself.'.format(aid, adjlist))

                if order[0] == '{':
                    order = order[1:-1].split(',')
                else:
                    order = [order]

                bonds[aid][aid2] = order

        # Check consistency using bonddict
        for atom1 in bonds:
            for atom2 in bonds[atom1]:
                if atom2 not in bonds:
                    raise InvalidAdjacencyListError(
                        'Error in adjacency list:\n{1}\nAtom {0:d} not in bond dictionary.'.format(atom2, adjlist))
                elif atom1 not in bonds[atom2]:
                    raise InvalidAdjacencyListError(
                        'Error in adjacency list:\n{2}\nFound bond between {0:d} and {1:d}, '
                        'but not the reverse'.format(atom1, atom2, adjlist))
                elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                    raise InvalidAdjacencyListError(
                        'Error in adjacency list: \n{4}\nFound bonds between {0:d} and {1:d}, but of different orders '
                        '"{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1], adjlist))

        # Convert bonddict to use Atom[group] and Bond[group] objects
        atomkeys = list(atomdict.keys())
        atomkeys.sort()
        for aid1 in atomkeys:
            atomkeys2 = list(bonds[aid1].keys())
            atomkeys2.sort()
            for aid2 in atomkeys2:
                if aid1 < aid2:
                    atom1 = atomdict[aid1]
                    atom2 = atomdict[aid2]
                    order = bonds[aid1][aid2]
                    if group:
                        bond = GroupBond(atom1, atom2, order)
                    elif len(order) == 1:
                        bond = Bond(atom1, atom2, order[0])
                    else:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nMultiple bond orders specified '
                                                        'for an atom.'.format(adjlist))
                    atom1.edges[atom2] = bond
                    atom2.edges[atom1] = bond

        if not group:
            if saturate_h:
                # Add explicit hydrogen atoms to complete structure if desired
                new_atoms = []
                for atom in atoms:
                    try:
                        valence = PeriodicSystem.valences[atom.symbol]
                    except KeyError:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nCannot add hydrogens: Unknown '
                                                        'valence for atom "{0}".'.format(atom.symbol, adjlist))
                    radical = atom.radical_electrons
                    order = atom.get_total_bond_order()
                    count = valence - radical - int(order) - 2 * (
                            atom.lone_pairs - PeriodicSystem.lone_pairs[atom.symbol])
                    for i in range(count):
                        a = Atom(element='H', radical_electrons=0, charge=0, label='', lone_pairs=0)
                        b = Bond(atom, a, 'S')
                        new_atoms.append(a)
                        atom.bonds[a] = b
                        a.bonds[atom] = b
                atoms.extend(new_atoms)

            # Calculate the multiplicity for the molecule and update the charges on each atom
            n_rad = 0  # total number of radical electrons
            for atom in atoms:
                atom.update_charge()
                n_rad += atom.radical_electrons
            multiplicity = n_rad + 1  # 2 s + 1, where s is the combined spin of unpaired electrons (s = 1/2 per unpaired electron)

        else:
            # Don't set a multiplicity for groups when converting from an old adjlist
            multiplicity = None

    except InvalidAdjacencyListError:
        logging.error("Troublesome adjacency list:\n" + adjlist)
        raise

    return atoms, multiplicity


###############################

re_intermediate_adjlist = re.compile(r'^\s*(\d*)\s+' +  # atom number digit
                                     r'(?P<label>\*\d*\s+)?' +  # optional label eg * or *2
                                     r'(?P<atomtype>\{?[A-Z]\S*)\s+' +  # atomtype eg R!H or {Cb,Cd}
                                     r'(?P<radicals>X|\d[STDQV]?|\{?\d[^}]*\})\s+' +  # radicals eg. X or 2T or {1,2,2T}
                                     r'(?P<lonepairs>\d)' +  # lone pairs eg. 0
                                     r'(?P<bonds>(\s+\{\d+\,(?:[SDTB]|\{.+?\})\},?)*)' +  # bonds, eg {2,S} {4,{S,D}}
                                     r'\s*$')  # the end!

re_old_adjlist = re.compile(r'^\s*(\d*)\s+' +  # atom number digit
                            r'(?P<label>\*\d*\s+)?' +  # optional label eg * or *2
                            r'(?P<atomtype>\{?[A-Z]\S*)\s+' +  # atomtype eg R!H or {Cb,Cd}
                            r'(?P<radicals>X|\d[STDQV]?|\{?\d[^}]*\})' +  # radicals eg. X or 2T or {1,2,2T}
                            r'(?P<bonds>(\s+\{\d+\,(?:[SDTB]|\{.+?\})\},?)*)' +  # bonds, eg {2,S} {4,{S,D}}
                            r'\s*$')  # the end!


def from_adjacency_list(adjlist, group=False, saturate_h=False):
    """
    Convert a string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects.
    """
    atoms = []
    atom_dict = {}
    bonds = {}
    multiplicity = None

    adjlist = adjlist.strip()
    lines = adjlist.splitlines()
    if adjlist == '' or len(lines) == 0:
        raise InvalidAdjacencyListError('Empty adjacency list.')

    # Detect old-style adjacency lists by looking at the last line's syntax
    last_line = lines[-1].strip()
    while not last_line:  # Remove any empty lines from the end
        lines.pop()
        last_line = lines[-1].strip()
    if re_intermediate_adjlist.match(last_line):
        logging.debug(
            "adjacency list:\n{1}\nline '{0}' looks like an intermediate style "
            "adjacency list".format(last_line, adjlist))
        return from_old_adjacency_list(adjlist, group=group, saturate_h=saturate_h)
    if re_old_adjlist.match(last_line):
        logging.debug(
            "Adjacency list:\n{1}\nline '{0}' looks like an old style adjacency list".format(last_line, adjlist))
        if not group:
            logging.debug("Will assume implicit H atoms")
        return from_old_adjacency_list(adjlist, group=group, saturate_h=(not group))

    # Interpret the first line if it contains a label
    if len(lines[0].split()) == 1:
        label = lines.pop(0)
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list.')

    # Interpret the second line if it contains a multiplicity
    if lines[0].split()[0] == 'multiplicity':
        line = lines.pop(0)
        if group:
            match = re.match(r'\s*multiplicity\s+\[\s*(\d(?:,\s*\d)*)\s*\]\s*$', line)
            if not match:
                rematch = re.match(r'\s*multiplicity\s+x\s*$', line)
                if not rematch:
                    raise InvalidAdjacencyListError("Invalid multiplicity line '{0}'. Should be a list like "
                                                    "'multiplicity [1,2,3]' or a wildcard 'multiplicity x'".format(line))
            else:
                # should match "multiplicity [1]" or " multiplicity   [ 1, 2, 3 ]" or " multiplicity [1,2,3]"
                # and whatever's inside the [] (excluding leading and trailing spaces) should be captured as group 1.
                # If a wildcard is desired, this line can be omitted or replaced with 'multiplicity x'
                # Multiplicities must be only one digit (i.e. less than 10)
                # The (?:,\s*\d)* matches patters like ", 2" 0 or more times, but doesn't capture them (because of the leading ?:)
                multiplicities = match.group(1).split(',')
                multiplicity = [int(i) for i in multiplicities]
        else:
            match = re.match(r'\s*multiplicity\s+\d+\s*$', line)
            if not match:
                raise InvalidAdjacencyListError("Invalid multiplicity line '{0}'. Should be an integer like "
                                                "'multiplicity 2'".format(line))
            multiplicity = int(line.split()[1])
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list: \n{0}'.format(adjlist))

    mistake1 = re.compile(r'\{[^}]*\s+[^}]*\}')
    # Iterate over the remaining lines, generating Atom or GroupAtom objects
    for line in lines:

        # Sometimes people put spaces after commas, which messes up the
        # parse-by-whitespace. Examples include '[Cd, Ct]'.
        if mistake1.search(line):
            raise InvalidAdjacencyListError(
                "{1} Shouldn't have spaces inside braces:\n{0}".format(mistake1.search(line).group(), adjlist)
            )

        # Sometimes commas are used to delimit bonds in the bond list,
        # so replace them just in case
        line = line.replace('},{', '} {')

        data = line.split()

        # Skip if blank line
        if len(data) == 0:
            continue

        # First item is index for atom
        # Sometimes these have a trailing period (as if in a numbered list),
        # so remove it just in case
        aid = int(data[0].strip('.'))

        # If second item starts with '*', then atom is labeled
        label = ''
        index = 1
        if data[1][0] == '*':
            label = data[1]
            index += 1

        # Next is the element or functional group element
        # A list can be specified with the {,} syntax
        atom_type = data[index]
        if atom_type[0] == '[':
            if not group:
                raise InvalidAdjacencyListError("Error on:\n{0}\nA molecule should not assign more than one "
                                                "atomtype per atom.".format(adjlist))
            atom_type = atom_type[1:-1].split(',')
        else:
            atom_type = [atom_type]
        index += 1

        # Next the number of unpaired electrons
        unpaired_electrons = []
        u_state = data[index]
        if u_state[0] == 'u':
            if u_state[1] == '[':
                u_state = u_state[2:-1].split(',')
            else:
                u_state = [u_state[1]]
            for u in u_state:
                if u == '0':
                    unpaired_electrons.append(0)
                elif u == '1':
                    unpaired_electrons.append(1)
                elif u == '2':
                    unpaired_electrons.append(2)
                elif u == '3':
                    unpaired_electrons.append(3)
                elif u == '4':
                    unpaired_electrons.append(4)
                elif u == 'x':
                    if not group:
                        raise InvalidAdjacencyListError("Error on:\n{0}\nA molecule should not assign a wildcard to "
                                                        "number of unpaired electrons.".format(adjlist))
                else:
                    raise InvalidAdjacencyListError('Number of unpaired electrons not recognized on\n{0}.'.format(adjlist))
            index += 1
        else:
            raise InvalidAdjacencyListError('Number of unpaired electrons not defined on\n{0}.'.format(adjlist))

        # Next the number of lone electron pairs (if provided)
        lone_pairs = []
        if len(data) > index:
            lp_state = data[index]
            if lp_state[0] == 'p':
                if lp_state[1] == '[':
                    lp_state = lp_state[2:-1].split(',')
                else:
                    lp_state = [lp_state[1]]
                for lp in lp_state:
                    if lp == '0':
                        lone_pairs.append(0)
                    elif lp == '1':
                        lone_pairs.append(1)
                    elif lp == '2':
                        lone_pairs.append(2)
                    elif lp == '3':
                        lone_pairs.append(3)
                    elif lp == '4':
                        lone_pairs.append(4)
                    elif lp == 'x':
                        if not group:
                            raise InvalidAdjacencyListError("Error in adjacency list:\n{0}\nA molecule should not have "
                                                            "a wildcard assigned to number of lone pairs.".format(adjlist))
                    else:
                        raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nNumber of lone electron pairs '
                                                        'not recognized.'.format(adjlist))
                index += 1
            else:
                if not group:
                    lone_pairs.append(0)
        else:
            if not group:
                lone_pairs.append(0)

        # Next the number of partial charges (if provided)
        partial_charges = []
        if len(data) > index:
            e_state = data[index]
            if e_state[0] == 'c':
                if e_state[1] == '[':
                    e_state = e_state[2:-1].split(',')
                else:
                    e_state = [e_state[1:]]
                for e in e_state:
                    if e == '0':
                        partial_charges.append(0)
                    elif e == '+1':
                        partial_charges.append(1)
                    elif e == '+2':
                        partial_charges.append(2)
                    elif e == '+3':
                        partial_charges.append(3)
                    elif e == '+4':
                        partial_charges.append(4)
                    elif e == '-1':
                        partial_charges.append(-1)
                    elif e == '-2':
                        partial_charges.append(-2)
                    elif e == '-3':
                        partial_charges.append(-3)
                    elif e == '-4':
                        partial_charges.append(-4)
                    elif e == 'x':
                        if not group:
                            raise InvalidAdjacencyListError("Error on adjacency list:\n{0}\nA molecule should not have "
                                                            "a wildcard assigned to number of charges.".format(adjlist))
                    else:
                        raise InvalidAdjacencyListError('Error on adjacency list:\n{0}\nNumber of partial charges '
                                                        'not recognized.'.format(adjlist))
                index += 1
            else:
                if not group:
                    partial_charges.append(0)
        else:
            if not group:
                partial_charges.append(0)

        # Next the isotope (if provided)
        isotope = -1
        if len(data) > index:
            i_state = data[index]
            if i_state[0] == 'i':
                isotope = int(i_state[1:])
                index += 1

        # Next ring membership info (if provided)
        props = {}
        if len(data) > index:
            r_state = data[index]
            if r_state[0] == 'r':
                props['inRing'] = bool(int(r_state[1]))
                index += 1

        # Create a new atom based on the above information
        if group:
            atom = GroupAtom(atom_type, unpaired_electrons, partial_charges, label, lone_pairs, props)
        else:
            atom = Atom(atom_type[0], unpaired_electrons[0], partial_charges[0], label, lone_pairs[0])
            if isotope != -1:
                atom.element = get_element(atom.number, isotope)

        # Add the atom to the list
        atoms.append(atom)
        atom_dict[aid] = atom

        # Process list of bonds
        bonds[aid] = {}
        for datum in data[index:]:

            # Sometimes commas are used to delimit bonds in the bond list,
            # so strip them just in case
            datum = datum.strip(',')

            aid2, comma, order = datum[1:-1].partition(',')
            aid2 = int(aid2)
            if aid == aid2:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAttempted to create a bond between '
                                                'atom {0:d} and itself.'.format(aid, adjlist))

            if order[0] == '[':
                order = order[1:-1].split(',')
            else:
                order = [order]

            bonds[aid][aid2] = order

    # Check consistency using bonddict
    for atom1 in bonds:
        for atom2 in bonds[atom1]:
            if atom2 not in bonds:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{1}\nAtom {0:d} not in bond '
                                                'dictionary.'.format(atom2, adjlist))
            elif atom1 not in bonds[atom2]:
                raise InvalidAdjacencyListError('Error in adjacency list:\n{2}\nFound bond between {0:d} and {1:d}, '
                                                'but not the reverse.'.format(atom1, atom2, adjlist))
            elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                raise InvalidAdjacencyListError(
                    'Error in adjacency list:\n{4}\nFound bonds between {0:d} and {1:d}, but of different orders '
                    '"{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1], adjlist))

    # Convert bonddict to use Atom[group] and Bond[group] objects
    atomkeys = list(atom_dict.keys())
    atomkeys.sort()
    for aid1 in atomkeys:
        atomkeys2 = list(bonds[aid1].keys())
        atomkeys2.sort()
        for aid2 in atomkeys2:
            if aid1 < aid2:
                atom1 = atom_dict[aid1]
                atom2 = atom_dict[aid2]
                order = bonds[aid1][aid2]
                if group:
                    bond = GroupBond(atom1, atom2, order)
                elif len(order) == 1:
                    bond = Bond(atom1, atom2, order[0])
                else:
                    raise InvalidAdjacencyListError('Error in adjacency list:\n{0}\nMultiple bond orders specified for '
                                                    'an atom in a Molecule.'.format(adjlist))
                atom1.edges[atom2] = bond
                atom2.edges[atom1] = bond

    if saturate_h:
        # Add explicit hydrogen atoms to complete structure if desired
        if not group:
            Saturator.saturate(atoms)

    # Consistency checks
    if not group:
        # Molecule consistency check
        # Electron and valency consistency check for each atom
        for atom in atoms:
            ConsistencyChecker.check_partial_charge(atom)

        n_rad = sum([atom.radical_electrons for atom in atoms])
        absolute_spin_per_electron = 1 / 2.
        if multiplicity is None:
            multiplicity = 2 * (n_rad * absolute_spin_per_electron) + 1

        ConsistencyChecker.check_multiplicity(n_rad, multiplicity)
        for atom in atoms:
            ConsistencyChecker.check_hund_rule(atom, multiplicity)
        return atoms, multiplicity
    else:
        # Currently no group consistency check
        return atoms, multiplicity


def to_adjacency_list(atoms, multiplicity, label=None, group=False, remove_h=False, remove_lone_pairs=False,
                      old_style=False):
    """
    Convert a chemical graph defined by a list of `atoms` into a string
    adjacency list.
    """
    if not atoms:
        return ''

    if old_style:
        return to_old_adjacency_list(atoms, multiplicity, label, group, remove_h)

    adjlist = ''

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if remove_h and all([atom.element.symbol == 'H' for atom in atoms]): remove_h = False
    except AttributeError:
        pass

    if label:
        adjlist += label + '\n'

    if group:
        if multiplicity:
            # Functional group should have a list of possible multiplicities.  
            # If the list is empty, then it does not need to be written
            adjlist += 'multiplicity [{0!s}]\n'.format(','.join(str(i) for i in multiplicity))
    else:
        assert isinstance(multiplicity, int), "Molecule should have an integer multiplicity"
        if multiplicity != 1 or any(atom.radical_electrons for atom in atoms):
            adjlist += 'multiplicity {0!r}\n'.format(multiplicity)

    # Determine the numbers to use for each atom
    atom_numbers = {}
    index = 0
    for atom in atoms:
        if remove_h and atom.element.symbol == 'H' and atom.label == '':
            continue
        atom_numbers[atom] = '{0:d}'.format(index + 1)
        index += 1

    atom_labels = dict([(atom, '{0}'.format(atom.label)) for atom in atom_numbers])

    atom_types = {}
    atom_unpaired_electrons = {}
    atom_lone_pairs = {}
    atom_charge = {}
    atom_isotope = {}
    atom_props = {}
    if group:
        for atom in atom_numbers:
            # Atom type(s)
            if len(atom.atomtype) == 1:
                atom_types[atom] = atom.atomtype[0].label
            else:
                atom_types[atom] = '[{0}]'.format(','.join([a.label for a in atom.atomtype]))
            # Unpaired Electron(s)
            if len(atom.radical_electrons) == 1:
                atom_unpaired_electrons[atom] = str(atom.radical_electrons[0])
            elif len(atom.radical_electrons) == 0:
                atom_unpaired_electrons[atom] = 'x'  # Empty list indicates wildcard
            else:
                atom_unpaired_electrons[atom] = '[{0}]'.format(','.join([str(radical) for radical in atom.radical_electrons]))

            # Lone Electron Pair(s)
            if len(atom.lone_pairs) == 1:
                atom_lone_pairs[atom] = str(atom.lone_pairs[0])
            elif len(atom.lone_pairs) == 0:
                atom_lone_pairs[atom] = None  # Empty list indicates wildcard
            else:
                atom_lone_pairs[atom] = '[{0}]'.format(','.join([str(pair) for pair in atom.lone_pairs]))

            # Charges
            if len(atom.charge) == 1:
                atom_charge[atom] = '+' + str(atom.charge[0]) if atom.charge[0] > 0 else str(atom.charge[0])
            elif len(atom.charge) == 0:
                atom_charge[atom] = None  # Empty list indicates wildcard
            else:
                atom_charge[atom] = '[{0}]'.format(','.join(['+'+str(charge) if charge > 0 else ''+str(charge) for charge in atom.charge]))

            # Isotopes
            atom_isotope[atom] = -1

            # Other props
            props = []
            if 'inRing' in atom.props:
                props.append(' r{0}'.format(int(atom.props['inRing'])))
            atom_props[atom] = props
    else:
        for atom in atom_numbers:
            # Atom type
            atom_types[atom] = '{0}'.format(atom.element.symbol)
            # Unpaired Electron(s)
            atom_unpaired_electrons[atom] = '{0}'.format(atom.radical_electrons)
            # Lone Electron Pair(s)
            atom_lone_pairs[atom] = str(atom.lone_pairs)
            # Partial Charge(s)
            atom_charge[atom] = '+' + str(atom.charge) if atom.charge > 0 else '' + str(atom.charge)
            # Isotopes
            atom_isotope[atom] = atom.element.isotope

    # Determine field widths
    atom_number_width = max([len(s) for s in atom_numbers.values()]) + 1
    atom_label_width = max([len(s) for s in atom_labels.values()])
    if atom_label_width > 0:
        atom_label_width += 1
    atom_type_width = max([len(s) for s in atom_types.values()]) + 1
    atom_unpaired_electrons_width = max([len(s) for s in atom_unpaired_electrons.values()])

    # Assemble the adjacency list
    for atom in atoms:
        if atom not in atom_numbers:
            continue

        # Atom number
        adjlist += '{0:<{1:d}}'.format(atom_numbers[atom], atom_number_width)
        # Atom label
        adjlist += '{0:<{1:d}}'.format(atom_labels[atom], atom_label_width)
        # Atom type(s)
        adjlist += '{0:<{1:d}}'.format(atom_types[atom], atom_type_width)
        # Unpaired Electron(s)
        adjlist += 'u{0:<{1:d}}'.format(atom_unpaired_electrons[atom], atom_unpaired_electrons_width)
        # Lone Electron Pair(s)
        if atom_lone_pairs[atom] is not None:
            adjlist += ' p{0}'.format(atom_lone_pairs[atom])
        # Partial charges
        if atom_charge[atom] is not None:
            adjlist += ' c{0}'.format(atom_charge[atom])
        # Isotopes
        if atom_isotope[atom] != -1:
            adjlist += ' i{0}'.format(atom_isotope[atom])
        if group and len(atom_props[atom]) > 0:
            for prop in atom_props[atom]:
                adjlist += prop

        # Bonds list
        atoms2 = list(atom.bonds.keys())
        # sort them the same way as the atoms
        atoms2.sort(key=atoms.index)

        for atom2 in atoms2:
            if atom2 not in atom_numbers:
                continue

            bond = atom.bonds[atom2]
            adjlist += ' {{{0},'.format(atom_numbers[atom2])

            # Bond type(s)
            if group:
                code = '[{0}]'
                if len(bond.order) == 1:
                    code = '{0}'
                # preference is for string representation, backs down to number
                # numbers if doesn't work
                try:
                    adjlist += code.format(','.join(bond.get_order_str()))
                except ValueError:
                    adjlist += code.format(','.join(str(bond.get_order_num())))
            else:
                # preference is for string representation, backs down to number
                # numbers if doesn't work
                try:
                    adjlist += bond.get_order_str()
                except ValueError:
                    adjlist += str(bond.get_order_num())
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist


def get_old_electron_state(atom):
    """
    Get the old adjacency list format electronic state
    """
    additional_lone_pairs = atom.lone_pairs - PeriodicSystem.lone_pairs[atom.element.symbol]
    electrons = atom.radical_electrons + additional_lone_pairs * 2
    if electrons == 0:
        electron_state = '0'
    elif electrons == 1:
        electron_state = '1'
    elif electrons == 2:
        if additional_lone_pairs == 0:
            electron_state = '2T'
        elif additional_lone_pairs == 1:
            electron_state = '2S'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    elif electrons == 3:
        if additional_lone_pairs == 0:
            electron_state = '3Q'
        elif additional_lone_pairs == 1:
            electron_state = '3D'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    elif electrons == 4:
        if additional_lone_pairs == 0:
            electron_state = '4V'
        elif additional_lone_pairs == 1:
            electron_state = '4T'
        elif additional_lone_pairs == 2:
            electron_state = '4S'
        else:
            raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    else:
        raise InvalidAdjacencyListError("Cannot find electron state of atom {0}".format(atom))
    return electron_state


def to_old_adjacency_list(atoms, multiplicity=None, label=None, group=False, remove_h=False):
    """
    Convert a chemical graph defined by a list of `atoms` into a string old-style 
    adjacency list that can be used in RMG-Java.  Currently not working for groups.
    """
    warnings.warn("The old adjacency lists are no longer supported and may be"
                  " removed in version 2.3.", DeprecationWarning)
    adjlist = ''

    if group:
        raise InvalidAdjacencyListError("Not yet implemented.")
    # Filter out all non-valid atoms
    if not group:
        for atom in atoms:
            if atom.element.symbol in ['He', 'Ne', 'Ar', 'N']:
                raise InvalidAdjacencyListError("Old-style adjacency list does not accept He, Ne, Ar, N elements.")

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if remove_h and all([atom.element.symbol == 'H' for atom in atoms]):
            remove_h = False
    except AttributeError:
        pass

    if label:
        adjlist += label + '\n'

    # Determine the numbers to use for each atom
    atom_numbers = {}
    index = 0
    for atom in atoms:
        if remove_h and atom.element.symbol == 'H' and atom.label == '': continue
        atom_numbers[atom] = '{0:d}'.format(index + 1)
        index += 1

    atom_labels = dict([(atom, '{0}'.format(atom.label)) for atom in atom_numbers])

    atom_types = {}
    atom_electron_states = {}
    if group:
        raise InvalidAdjacencyListError("Not yet implemented.")
    else:
        for atom in atom_numbers:
            # Atom type
            atom_types[atom] = '{0}'.format(atom.element.symbol)
            # Electron state(s)
            atom_electron_states[atom] = '{0}'.format(get_old_electron_state(atom))

    # Determine field widths
    atom_number_width = max([len(s) for s in atom_numbers.values()]) + 1
    atom_label_width = max([len(s) for s in atom_labels.values()])
    if atom_label_width > 0:
        atom_label_width += 1
    atom_type_width = max([len(s) for s in atom_types.values()]) + 1
    atom_electron_state_width = max([len(s) for s in atom_electron_states.values()])

    # Assemble the adjacency list
    for atom in atoms:
        if atom not in atom_numbers:
            continue

        # Atom number
        adjlist += '{0:<{1:d}}'.format(atom_numbers[atom], atom_number_width)
        # Atom label
        adjlist += '{0:<{1:d}}'.format(atom_labels[atom], atom_label_width)
        # Atom type(s)
        adjlist += '{0:<{1:d}}'.format(atom_types[atom], atom_type_width)
        # Electron state(s)
        adjlist += '{0:<{1:d}}'.format(atom_electron_states[atom], atom_electron_state_width)

        # Bonds list
        atoms2 = list(atom.bonds.keys())
        # sort them the same way as the atoms
        atoms2.sort(key=atoms.index)

        for atom2 in atoms2:
            if atom2 not in atom_numbers:
                continue

            bond = atom.bonds[atom2]
            adjlist += ' {{{0},'.format(atom_numbers[atom2])

            # Bond type(s)
            if group:
                if len(bond.order) == 1:
                    adjlist += bond.get_order_str()[0]
                else:
                    adjlist += '{{{0}}}'.format(','.join(bond.get_order_str()))
            else:
                adjlist += bond.get_order_str()
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist
