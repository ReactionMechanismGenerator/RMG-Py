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
This module provides classes and methods for working with molecules and
molecular configurations. A molecule is represented internally using a graph
data type, where atoms correspond to vertices and bonds correspond to edges.
Both :class:`Atom` and :class:`Bond` objects store semantic information that
describe the corresponding atom or bond.
"""

import itertools
import logging
import os
from collections import OrderedDict, defaultdict
from copy import deepcopy
from urllib.parse import quote
from operator import attrgetter

import cython
import numpy as np

import rmgpy.constants as constants
import rmgpy.molecule.converter as converter
import rmgpy.molecule.element as elements
import rmgpy.molecule.group as gr
import rmgpy.molecule.resonance as resonance
import rmgpy.molecule.translator as translator
from rmgpy.exceptions import DependencyError
from rmgpy.molecule.adjlist import Saturator
from rmgpy.molecule.atomtype import AtomType, ATOMTYPES, get_atomtype, AtomTypeError
from rmgpy.molecule.element import bdes
from rmgpy.molecule.graph import Vertex, Edge, Graph, get_vertex_connectivity_value
from rmgpy.molecule.kekulize import kekulize
from rmgpy.molecule.pathfinder import find_shortest_path
from rmgpy.molecule.fragment import CuttingLabel

################################################################################

# helper function for sorting
def _skip_first(in_tuple):
    return in_tuple[1:]

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
    `site`               ``str``             type of adsorption site
    `morphology`         ``str``             morphology of the adsorption site
    ==================== =================== ====================================

    Additionally, the ``mass``, ``number``, and ``symbol`` attributes of the
    atom's element can be read (but not written) directly from the atom object,
    e.g. ``atom.symbol`` instead of ``atom.element.symbol``.
    """

    def __init__(self, element=None, radical_electrons=0, charge=0, label='', lone_pairs=-100, site='', morphology='', 
                 coords=np.array([]), id=-1, props=None):
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
        self.site = site 
        self.morphology = morphology
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
            'site': self.site,
            'morphology': self.morphology,
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
        self.site = d['site']
        self.morphology = d['morphology']

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
        if issubclass(type(other), Vertex):
            return self.sorting_key < other.sorting_key
        else:
            raise NotImplementedError('Cannot perform less than comparison between Atom and '
                                      '{0}.'.format(type(other).__name__))

    def __gt__(self, other):
        """Define greater than comparison. For comparing against other Atom objects (e.g. when sorting)."""
        if issubclass(type(other), Vertex):
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
        cython.declare(atom=Atom, ap=gr.GroupAtom)
        if isinstance(other, Atom):
            atom = other
            if strict:
                return (self.element is atom.element
                        and self.radical_electrons == atom.radical_electrons
                        and self.lone_pairs == atom.lone_pairs
                        and self.charge == atom.charge
                        and self.atomtype is atom.atomtype
                        and self.site == atom.site
                        and self.morphology == atom.morphology)

            else:
                return self.element is atom.element
        elif isinstance(other, gr.GroupAtom):
            cython.declare(a=AtomType, radical=cython.short, lp=cython.short, charge=cython.short)
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
            if ap.site:
                for site in ap.site:
                    if self.site == site: break
                else:
                    return False
            if ap.morphology:
                for morphology in ap.morphology:
                    if self.morphology == morphology: break
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
            cython.declare(atom=gr.GroupAtom, a=AtomType, radical=cython.short, lp=cython.short, charge=cython.short)
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
            if atom.site:
                for site in atom.site:
                    if self.site == site: break
                else:
                    return False
            if atom.morphology:
                for morphology in atom.morphology:
                    if self.morphology == morphology: break
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
        cython.declare(a=Atom)
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
        a.site = self.site
        a.morphology = self.morphology
        a.coords = self.coords[:]
        a.id = self.id
        a.props = deepcopy(self.props)
        return a

    def is_electron(self):
        """
        Return ``True`` if the atom represents an electron or ``False`` if
        not.
        """
        return self.element.number == -1
    
    def is_proton(self):
        """
        Return ``True`` if the atom represents a proton or ``False`` if
        not.
        """

        if self.element.number == 1 and self.charge == 1:
            return True
        return False

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

    def is_lithium(self):
        """
        Return ``True`` if the atom represents a hydrogen atom or ``False`` if
        not.
        """
        return self.element.number == 3

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

    def is_bonded_to_surface(self):
        """
        Return ``True`` if the atom is bonded to a surface atom `X`
        ``False`` if it is not
        """
        cython.declare(bonded_atom=Atom)
        for bonded_atom in self.bonds.keys():
            if bonded_atom.is_surface_site():
                return True
        return False

    def is_bonded_to_halogen(self):
        """
        Return ``True`` if the atom is bonded to at least one halogen (F, Cl, Br, or I)
        ``False`` if it is not
        """
        cython.declare(bonded_atom=Atom)
        for bonded_atom in self.bonds.keys():
            if bonded_atom.is_halogen():
                return True
        return False

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
        cython.declare(radical_electrons=cython.short)
        # Set the new radical electron count
        radical_electrons = self.radical_electrons = self.radical_electrons - 1
        if radical_electrons < 0:
            raise gr.ActionError('Unable to update Atom due to LOSE_RADICAL action: '
                                 'Invalid radical electron set "{0}".'.format(self.radical_electrons))

    def increment_charge(self):
        """
        Update the atom pattern as a result of applying a GAIN_CHARGE action
        """
        self.charge += 1

    def decrement_charge(self):
        """
        Update the atom pattern as a result of applying a LOSE_CHARGE action
        """
        self.charge -= 1

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
        if self.is_electron():
            self.charge = -1
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
        elif act == 'GAIN_CHARGE':
            for i in range(action[2]): self.increment_charge()
        elif act == 'LOSE_CHARGE':
            for i in range(abs(action[2])): self.decrement_charge()
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
        cython.declare(bond=Bond, bp=gr.GroupBond)
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
        elif self.is_reaction_bond():
            return 'R'
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
        elif new_order == 'R':
            self.order = 0.05
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
        cython.declare(b=Bond)
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

    def is_double_or_triple(self):
        """
        Return ``True`` if the bond represents a double or triple bond or ``False``
        if not.
        """
        return self.is_order(2) or self.is_order(3)

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

    def is_reaction_bond(self):
        """
        Return ``True`` if the bond represents a reaction bond or ``False`` if
        not.
        """
        return self.is_order(0.05)
    
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
        bond_symbol_mapping = {0.05: '~', 0.1: '~', 1: '-', 1.5: ':', 2: '=', 3: '#'}
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


class Molecule(Graph):
    """
    A representation of a molecular structure using a graph data type, extending
    the :class:`Graph` class. Attributes are:

    ======================= =========== ========================================
    Attribute               Type        Description
    ======================= =========== ========================================
    `atoms`                 ``list``    A list of Atom objects in the molecule
    `symmetry_number`       ``float``   The (estimated) external + internal symmetry number of the molecule, modified for chirality
    `multiplicity`          ``int``     The multiplicity of this species, multiplicity = 2*total_spin+1
    `reactive`              ``bool``    ``True`` (by default) if the molecule participates in reaction families.
                                            It is set to ``False`` by the filtration functions if a non
                                            representative resonance structure was generated by a template reaction
    `props`                 ``dict``    A list of properties describing the state of the molecule.
    `inchi`                 ``str``     A string representation of the molecule in InChI
    `smiles`                ``str``     A string representation of the molecule in SMILES
    `fingerprint`           ``str``     A representation for fast comparison, set as molecular formula
    `metal`                 ``str``     The metal of the metal surface the molecule is associated with
    `facet`                 ``str``     The facet of the metal surface the molecule is associated with
    ======================= =========== ========================================

    A new molecule object can be easily instantiated by passing the `smiles` or
    `inchi` string representing the molecular structure.
    """

    def __init__(self, atoms=None, symmetry=-1, multiplicity=-187, reactive=True, props=None, inchi='', smiles='', 
                 metal='', facet=''):
        Graph.__init__(self, atoms)
        self.symmetry_number = symmetry
        self.multiplicity = multiplicity
        self.reactive = reactive
        self._fingerprint = None
        self._inchi = None
        self._smiles = None
        self.props = props or {}
        self.metal = metal
        self.facet = facet

        if inchi and smiles:
            logging.warning('Both InChI and SMILES provided for Molecule instantiation, '
                            'using InChI and ignoring SMILES.')
        if inchi:
            self.from_inchi(inchi)
            self._inchi = inchi
        elif smiles:
            self.from_smiles(smiles)
            self._smiles = smiles

        if multiplicity != -187:  # it was set explicitly, so re-set it (from_smiles etc may have changed it)
            self.multiplicity = multiplicity

    def __deepcopy__(self, memo):
        return self.copy(deep=True)

    def __hash__(self):
        """
        Define a custom hash method to allow Molecule objects to be used in dictionaries and sets.

        Use the fingerprint property, which is currently defined as the molecular formula, though
        this is not an ideal hash since there will be significant hash collision, leading to inefficient lookups.
        """
        return hash(('Molecule', self.fingerprint))

    def __eq__(self, other):
        """Method to test equality of two Molecule objects."""
        return self is other or (isinstance(other, Molecule) and
                                 self.fingerprint == other.fingerprint and
                                 self.is_isomorphic(other))

    def __lt__(self, other):
        """Define less than comparison. For comparing against other Molecule objects (e.g. when sorting)."""
        if isinstance(other, Molecule):
            return self.sorting_key < other.sorting_key
        else:
            raise NotImplementedError('Cannot perform less than comparison between Molecule and '
                                      '{0}.'.format(type(other).__name__))

    def __gt__(self, other):
        """Define greater than comparison. For comparing against other Molecule objects (e.g. when sorting)."""
        if isinstance(other, Molecule):
            return self.sorting_key > other.sorting_key
        else:
            raise NotImplementedError('Cannot perform greater than comparison between Molecule and '
                                      '{0}.'.format(type(other).__name__))

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return '<Molecule "{0}">'.format(self.to_smiles())

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        cython.declare(multiplicity=cython.int)
        multiplicity = self.multiplicity
        try:
            if multiplicity != self.get_radical_count() + 1:
                return 'Molecule(smiles="{0}", multiplicity={1:d})'.format(self.to_smiles(), multiplicity)
            return 'Molecule(smiles="{0}")'.format(self.to_smiles())
        except KeyError:
            logging.warning('Could not generate SMILES for this molecule object.'
                            ' Likely due to a keyerror when converting to RDKit'
                            ' Here is molecules AdjList: {}'.format(self.to_adjacency_list()))
            return 'Molecule().from_adjacency_list"""{}"""'.format(self.to_adjacency_list())

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Molecule, (self.vertices, self.symmetry_number, self.multiplicity, self.reactive, self.props, self.metal, self.facet))

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

    @property
    def fingerprint(self):
        """
        Fingerprint used to accelerate graph isomorphism comparisons with
        other molecules. The fingerprint is a short string containing a
        summary of selected information about the molecule. Two fingerprint
        strings matching is a necessary (but not sufficient) condition for
        the associated molecules to be isomorphic.

        Use an expanded molecular formula to also enable sorting.
        """
        if self._fingerprint is None:
            # Include these elements in this order at minimum
            element_dict = {'C': 0, 'H': 0, 'N': 0, 'O': 0, 'S': 0}
            all_elements = sorted(self.get_element_count().items(), key=lambda x: x[0])  # Sort alphabetically
            element_dict.update(all_elements)
            self._fingerprint = ''.join([f'{symbol}{num:0>2}' for symbol, num in element_dict.items()])
        return self._fingerprint

    @fingerprint.setter
    def fingerprint(self, fingerprint):
        self._fingerprint = fingerprint

    @property
    def inchi(self):
        """InChI string for this molecule. Read-only."""
        if self._inchi is None:
            self._inchi = self.to_inchi()
        return self._inchi

    @property
    def smiles(self):
        """SMILES string for this molecule. Read-only."""
        if self._smiles is None:
            self._smiles = self.to_smiles()
        return self._smiles

    @property
    def sorting_key(self):
        """Returns a sorting key for comparing Molecule objects. Read-only"""
        return self.fingerprint

    def add_atom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.add_vertex(atom)

    def add_bond(self, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.add_edge(bond)

    def get_bonds(self, atom):
        """
        Return a dictionary of the bonds involving the specified `atom`.
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
        Returns ``True`` iff the molecule contains an 'X' surface site.
        """
        cython.declare(atom=Atom)
        for atom in self.atoms:
            if atom.symbol == 'X':
                return True
        return False

    def number_of_surface_sites(self):
        """
        Returns the number of surface sites in the molecule.
        e.g. 2 for a bidentate adsorbate
        """
        cython.declare(atom=Atom)
        cython.declare(count=cython.int)
        count = 0
        for atom in self.atoms:
            if atom.is_surface_site():
                count += 1
        return count

    def is_surface_site(self):
        """Returns ``True`` iff the molecule is nothing but a surface site 'X'."""
        return len(self.atoms) == 1 and self.atoms[0].is_surface_site()

    def is_electron(self):
        """Returns ``True`` iff the molecule is nothing but an electron 'e'."""
        return len(self.atoms) == 1 and self.atoms[0].is_electron()

    def is_proton(self):
        """Returns ``True`` iff the molecule is nothing but a proton 'H+'."""
        return len(self.atoms) == 1 and self.atoms[0].is_proton()

    def remove_atom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.remove_vertex(atom)

    def remove_bond(self, bond):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        self._fingerprint = self._inchi = self._smiles = None
        return self.remove_edge(bond)

    def remove_van_der_waals_bonds(self):
        """
        Remove all van der Waals bonds.
        """
        cython.declare(bond=Bond)
        for bond in self.get_all_edges():
            if bond.is_van_der_waals():
                self.remove_bond(bond)

    def sort_atoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.

        This function orders atoms using several attributes in atom.getDescriptor().
        Currently it sorts by placing heaviest atoms first and hydrogen atoms last.
        Placing hydrogens last during sorting ensures that functions with hydrogen
        removal work properly.
        """
        cython.declare(vertex=Vertex, a=Atom, index=int)
        for vertex in self.vertices:
            if vertex.sorting_label < 0:
                self.update_connectivity_values()
                break
        self.atoms.sort(reverse=True)
        for index, vertex in enumerate(self.vertices):
            vertex.sorting_label = index

    def update_charge(self):

        for atom in self.atoms:
            if not isinstance(atom, CuttingLabel):
                atom.update_charge()

    def update(self, log_species=True, raise_atomtype_exception=True, sort_atoms=True):
        """
        Update the lone_pairs, charge, and atom types of atoms.
        Update multiplicity, and sort atoms (if ``sort_atoms`` is ``True``)
        Does not necessarily update the connectivity values (which are used in isomorphism checks)
        If you need that, call update_connectivity_values()
        """

        self.update_lone_pairs()
        self.update_charge()
        self.update_atomtypes(log_species=log_species, raise_exception=raise_atomtype_exception)
        self.update_multiplicity()
        if sort_atoms:
            self.sort_atoms()
        self.identify_ring_membership()

    def get_formula(self):
        """
        Return the molecular formula for the molecule.
        """
        cython.declare(atom=Atom, symbol=str, elements=dict, keys=list, formula=str)
        cython.declare(hasCarbon=cython.bint, hasHydrogen=cython.bint)

        # Count the number of each element in the molecule
        element_dict = {}
        for atom in self.vertices:
            symbol = atom.element.symbol
            element_dict[symbol] = element_dict.get(symbol, 0) + 1

        # Use the Hill system to generate the formula.
        # If you change this algorithm consider also updating 
        # the chemkin.write_thermo_entry method
        formula = ''

        # Carbon and hydrogen always come first if carbon is present
        if 'C' in element_dict:
            count = element_dict['C']
            formula += 'C{0:d}'.format(count) if count > 1 else 'C'
            del element_dict['C']
            if 'H' in element_dict:
                count = element_dict['H']
                formula += 'H{0:d}'.format(count) if count > 1 else 'H'
                del element_dict['H']

        # Other atoms are in alphabetical order
        # (This includes hydrogen if carbon is not present)
        keys = sorted(element_dict.keys())
        for key in keys:
            count = element_dict[key]
            formula += '{0}{1:d}'.format(key, count) if count > 1 else key

        return formula

    def get_molecular_weight(self):
        """
        Return the molecular weight of the molecule in kg/mol.
        """
        cython.declare(atom=Atom, mass=cython.double)
        mass = 0
        for atom in self.vertices:
            mass += atom.element.mass
        return mass

    def get_radical_count(self):
        """
        Return the total number of radical electrons on all atoms in the
        molecule. In this function, monoradical atoms count as one, biradicals
        count as two, etc.
        """
        cython.declare(atom=Atom, radicals=cython.short)
        radicals = 0
        for atom in self.vertices:
            radicals += atom.radical_electrons
        return radicals

    def get_singlet_carbene_count(self):
        """
        Return the total number of singlet carbenes (lone pair on a carbon atom)
        in the molecule. Counts the number of carbon atoms with a lone pair.
        In the case of [C] with two lone pairs, this method will return 1.
        """
        cython.declare(atom=Atom, carbenes=cython.short)
        carbenes = 0
        for atom in self.vertices:
            if atom.is_carbon() and atom.lone_pairs > 0:
                carbenes += 1
        return carbenes

    def get_num_atoms(self, element=None):
        """
        Return the number of atoms in molecule.  If element is given, ie. "H" or "C",
        the number of atoms of that element is returned.
        """
        cython.declare(numAtoms=cython.int, atom=Atom)
        if element is None:
            return len(self.vertices)
        else:
            num_atoms = 0
            for atom in self.vertices:
                if atom.element.symbol == element:
                    num_atoms += 1
            return num_atoms

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        cython.declare(g=Graph, other=Molecule, i=int, v1=Vertex, v2=Vertex)
        g = Graph.copy(self, deep)
        other = Molecule(g.vertices)
        # Copy connectivity values and sorting labels
        for i in range(len(self.vertices)):
            v1 = self.vertices[i]
            v2 = other.vertices[i]
            v2.connectivity1 = v1.connectivity1
            v2.connectivity2 = v1.connectivity2
            v2.connectivity3 = v1.connectivity3
            v2.sorting_label = v1.sorting_label
        other.multiplicity = self.multiplicity
        other.reactive = self.reactive
        other.metal = self.metal 
        other.facet = self.facet
        return other

    def merge(self, other):
        """
        Merge two molecules so as to store them in a single :class:`Molecule`
        object. The merged :class:`Molecule` object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Molecule(atoms=g.vertices)
        return molecule

    def split(self):
        """
        Convert a single :class:`Molecule` object containing two or more
        unconnected molecules into separate class:`Molecule` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Molecule(atoms=g.vertices)
            molecules.append(molecule)
        return molecules

    def delete_hydrogens(self):
        """
        Irreversibly delete all non-labeled hydrogens without updating
        connectivity values. If there's nothing but hydrogens, it does nothing.
        It destroys information; be careful with it.
        """
        cython.declare(atom=Atom, hydrogens=list)
        # Check that the structure contains at least one heavy atom
        for atom in self.vertices:
            if not atom.is_hydrogen():
                break
        else:
            # No heavy atoms, so leave explicit
            return
        hydrogens = []
        for atom in self.vertices:
            if atom.is_hydrogen() and atom.label == '':
                hydrogens.append(atom)
        # Remove the hydrogen atoms from the structure
        for atom in hydrogens:
            self.remove_atom(atom)

    def connect_the_dots(self, critical_distance_factor=0.45, raise_atomtype_exception=True):
        """
        Delete all bonds, and set them again based on the Atoms' coords.
        Does not detect bond type.
        """
        cython.declare(criticalDistance=float, i=int, atom1=Atom, atom2=Atom,
                       bond=Bond, atoms=list, zBoundary=float)
        # groupBond=GroupBond,
        self._fingerprint = None

        atoms = self.vertices

        # Ensure there are coordinates to work with
        for atom in atoms:
            assert len(atom.coords) != 0

        # If there are any bonds, remove them
        for atom1 in atoms:
            for bond in self.get_bonds(atom1):
                self.remove_edge(bond)

        # Sort atoms by distance on the z-axis
        sorted_atoms = sorted(atoms, key=lambda x: x.coords[2])

        for i, atom1 in enumerate(sorted_atoms):
            for atom2 in sorted_atoms[i + 1:]:
                # Set upper limit for bond distance
                critical_distance = (
                    atom1.element.cov_radius + atom2.element.cov_radius + critical_distance_factor) ** 2

                # First atom that is more than 4.0 Anstroms away in the z-axis, break the loop
                # Atoms are sorted along the z-axis, so all following atoms should be even further
                z_boundary = (atom1.coords[2] - atom2.coords[2]) ** 2
                if z_boundary > 16.0:
                    break

                distance_squared = sum((atom1.coords - atom2.coords) ** 2)

                if distance_squared > critical_distance or distance_squared < 0.40:
                    continue
                else:
                    # groupBond = GroupBond(atom1, atom2, [1,2,3,1.5])
                    bond = Bond(atom1, atom2, 1)
                    self.add_bond(bond)
        self.update_atomtypes(raise_exception=raise_atomtype_exception)

    def update_atomtypes(self, log_species=True, raise_exception=True):
        """
        Iterate through the atoms in the structure, checking their atom types
        to ensure they are correct (i.e. accurately describe their local bond
        environment) and complete (i.e. are as detailed as possible).
        
        If `raise_exception` is `False`, then the generic atomtype 'R' will
        be prescribed to any atom when get_atomtype fails. Currently used for
        resonance hybrid atom types.
        """
        # Because we use lonepairs to match atomtypes and default is -100 when unspecified,
        # we should update before getting the atomtype.
        self.update_lone_pairs()

        for atom in self.vertices:
            try:
                atom.atomtype = get_atomtype(atom, atom.edges)
            except AtomTypeError:
                if log_species:
                    logging.error("Could not update atomtypes for this molecule:\n{0}".format(self.to_adjacency_list()))
                if raise_exception:
                    raise
                atom.atomtype = ATOMTYPES['R']

    def update_multiplicity(self):
        """
        Update the multiplicity of a newly formed molecule.
        """
        # Assume this is always true
        # There are cases where 2 radical_electrons is a singlet, but
        # the triplet is often more stable, 
        self.multiplicity = self.get_radical_count() + 1

    def clear_labeled_atoms(self):
        """
        Remove the labels from all atoms in the molecule.
        """
        for atom in self.vertices:
            atom.label = ''

    def contains_labeled_atom(self, label):
        """
        Return :data:`True` if the molecule contains an atom with the label
        `label` and :data:`False` otherwise.
        """
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def get_labeled_atoms(self, label):
        """
        Return the atoms in the molecule that are labeled.
        """
        alist = [atom for atom in self.vertices if atom.label == label]
        if alist == []:
            raise ValueError(
                'No atom in the molecule \n{1}\n has the label "{0}".'.format(label, self.to_adjacency_list()))
        return alist

    def get_all_labeled_atoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
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
        """
        cython.declare(atom=Atom, element_count=dict, symbol=str, key=str)
        element_count = {}
        for atom in self.atoms:
            symbol = atom.element.symbol
            if symbol in element_count:
                element_count[symbol] += 1
            else:
                element_count[symbol] = 1

        return element_count

    def is_isomorphic(self, other, initial_map=None, generate_initial_map=False, save_order=False, strict=True):
        """
        Returns :data:`True` if two graphs are isomorphic and :data:`False`
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Molecule` object, or a :class:`TypeError` is raised.
        Also ensures multiplicities are also equal.

        Args:
            initial_map (dict, optional):          initial atom mapping to use
            generate_initial_map (bool, optional): if ``True``, initialize map by pairing atoms with same labels
            save_order (bool, optional):           if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):               if ``False``, perform isomorphism ignoring electrons
        """
        # It only makes sense to compare a Molecule to a Molecule for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Graph):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the quick isomorphism comparison using the fingerprint
        # Two fingerprint strings matching is a necessary (but not
        # sufficient!) condition for the associated molecules to be isomorphic
        if self.fingerprint != other.fingerprint:
            return False
        # check multiplicity
        if self.multiplicity != other.multiplicity:
            return False
        #check metal
        if self.metal != other.metal:
            return False 
        #check facet
        if self.facet != other.facet:
            return False
        # if given an initial map, ensure that it's valid.
        if initial_map:
            if not self.is_mapping_valid(other, initial_map, equivalent=True):
                return False

        # Do the full isomorphism comparison
        result = Graph.is_isomorphic(self, other, initial_map, generate_initial_map, save_order=save_order, strict=strict)
        return result

    def find_isomorphism(self, other, initial_map=None, save_order=False, strict=True):
        """
        Returns :data:`True` if `other` is isomorphic and :data:`False`
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`Molecule` object, or a :class:`TypeError` is raised.

        Args:
            initial_map (dict, optional): initial atom mapping to use
            save_order (bool, optional):  if ``True``, reset atom order after performing atom isomorphism
            strict (bool, optional):      if ``False``, perform isomorphism ignoring electrons
        """
        # It only makes sense to compare a Molecule to a Molecule for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Molecule):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))
        # Do the quick isomorphism comparison using the fingerprint
        # Two fingerprint strings matching is a necessary (but not
        # sufficient!) condition for the associated molecules to be isomorphic
        if self.fingerprint != other.fingerprint:
            return []
        # check multiplicity
        if self.multiplicity != other.multiplicity:
            return []
        #check metal
        if self.metal != other.metal:
            return []
        #check facet
        if self.facet != other.facet:
            return []
        # Do the isomorphism comparison
        result = Graph.find_isomorphism(self, other, initial_map, save_order=save_order, strict=strict)
        return result

    def is_subgraph_isomorphic(self, other, initial_map=None, generate_initial_map=False, save_order=False):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. The `initial_map` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        cython.declare(group=gr.Group, atom=Atom)
        cython.declare(carbonCount=cython.short, nitrogenCount=cython.short, oxygenCount=cython.short,
                       sulfurCount=cython.short, radicalCount=cython.short)
        cython.declare(L=list)
        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, gr.Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other

        # Check multiplicity
        if group.multiplicity:
            if self.multiplicity not in group.multiplicity: return False
        #check metal
        if group.metal:
            if self.metal not in group.metal: return False
        #check facet
        if group.facet:
            if self.facet not in group.facet: return False
        # Compare radical counts
        if self.get_radical_count() < group.radicalCount:
            return False

        # Compare element counts
        element_count = self.get_element_count()
        for element, count in group.elementCount.items():
            if element not in element_count:
                return False
            elif element_count[element] < count:
                return False

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
                    if len(set(atmlist)) != len(atmlist):  # skip entries that map multiple graph atoms to the same subgraph atom
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

        # Do the isomorphism comparison
        result = Graph.is_subgraph_isomorphic(self, other, initial_map, save_order=save_order)
        return result

    def find_subgraph_isomorphisms(self, other, initial_map=None, save_order=False):
        """
        Returns :data:`True` if `other` is subgraph isomorphic and :data:`False`
        otherwise. Also returns the lists all of valid mappings. The
        `initial_map` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        cython.declare(group=gr.Group, atom=Atom)
        cython.declare(carbonCount=cython.short, nitrogenCount=cython.short, oxygenCount=cython.short,
                       sulfurCount=cython.short, radicalCount=cython.short)

        # It only makes sense to compare a Molecule to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, gr.Group):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        group = other

        # Check multiplicity
        if group.multiplicity:
            if self.multiplicity not in group.multiplicity: return []
        #check metal
        if group.metal:
            if self.metal not in group.metal: return False
        #check facet
        if group.facet:
            if self.facet not in group.facet: return False
        # Compare radical counts
        if self.get_radical_count() < group.radicalCount:
            return []

        # Compare element counts
        element_count = self.get_element_count()
        for element, count in group.elementCount.items():
            if element not in element_count:
                return []
            elif element_count[element] < count:
                return []

        # Do the isomorphism comparison
        result = Graph.find_subgraph_isomorphisms(self, other, initial_map, save_order=save_order)
        return result

    def is_atom_in_cycle(self, atom):
        """
        Return :data:``True`` if ``atom`` is in one or more cycles in the structure,
        and :data:``False`` if not.
        """
        return self.is_vertex_in_cycle(atom)

    def is_bond_in_cycle(self, bond):
        """
        Return :data:``True`` if the bond between atoms ``atom1`` and ``atom2``
        is in one or more cycles in the graph, or :data:``False`` if not.
        """
        return self.is_edge_in_cycle(bond)

    def draw(self, path):
        """
        Generate a pictorial representation of the chemical graph using the
        :mod:`draw` module. Use `path` to specify the file to save
        the generated image to; the image type is automatically determined by
        extension. Valid extensions are ``.png``, ``.svg``, ``.pdf``, and
        ``.ps``; of these, the first is a raster format and the remainder are
        vector formats.
        """
        from .draw import MoleculeDrawer
        img_format = os.path.splitext(path)[-1][1:].lower()
        MoleculeDrawer().draw(self, img_format, target=path)

    def _repr_png_(self):
        """
        Return a png picture of the molecule, useful for ipython-qtconsole.
        """
        from .draw import MoleculeDrawer
        temp_file_name = 'temp_molecule.png'
        MoleculeDrawer().draw(self, 'png', temp_file_name)
        png = open(temp_file_name, 'rb').read()
        os.unlink(temp_file_name)
        return png

    def from_inchi(self, inchistr, backend='openbabel-first', raise_atomtype_exception=True):
        """
        Convert an InChI string `inchistr` to a molecular structure.

        RDKit and Open Babel are the two backends used in RMG. It is possible to use a
        single backend or try different backends in sequence. The available options for the ``backend``
        argument: 'openbabel-first'(default), 'rdkit-first', 'rdkit', or 'openbabel'.
        """
        translator.from_inchi(self, inchistr, backend, raise_atomtype_exception=raise_atomtype_exception)
        return self

    def from_augmented_inchi(self, aug_inchi, raise_atomtype_exception=True):
        """
        Convert an Augmented InChI string `aug_inchi` to a molecular structure.
        """
        translator.from_augmented_inchi(self, aug_inchi, raise_atomtype_exception=raise_atomtype_exception)
        return self

    def from_smiles(self, smilesstr, backend='openbabel-first', raise_atomtype_exception=True):
        """
        Convert a SMILES string `smilesstr` to a molecular structure.

        RDKit and Open Babel are the two backends used in RMG. It is possible to use a
        single backend or try different backends in sequence. The available options for the ``backend``
        argument: 'openbabel-first'(default), 'rdkit-first', 'rdkit', or 'openbabel'.
        """
        translator.from_smiles(self, smilesstr, backend, raise_atomtype_exception=raise_atomtype_exception)
        return self

    def from_smarts(self, smartsstr, raise_atomtype_exception=True):
        """
        Convert a SMARTS string `smartsstr` to a molecular structure. Uses
        `RDKit <https://rdkit.org/>`_ to perform the conversion.
        This Kekulizes everything, removing all aromatic atom types.
        """
        translator.from_smarts(self, smartsstr, raise_atomtype_exception=raise_atomtype_exception)
        return self

    def from_adjacency_list(self, adjlist, saturate_h=False, raise_atomtype_exception=True,
                            raise_charge_exception=False, check_consistency=True):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        from rmgpy.molecule.adjlist import from_adjacency_list

        self.vertices, self.multiplicity, self.metal, self.facet = from_adjacency_list(adjlist, group=False, saturate_h=saturate_h,
                                                               check_consistency=check_consistency)
        self.update_atomtypes(raise_exception=raise_atomtype_exception)
        self.identify_ring_membership()

        # Check if multiplicity is possible
        n_rad = self.get_radical_count()
        multiplicity = self.multiplicity
        if not (n_rad + 1 == multiplicity or n_rad - 1 == multiplicity or
                n_rad - 3 == multiplicity or n_rad - 5 == multiplicity):
            raise ValueError('Impossible multiplicity for molecule\n{0}\n multiplicity = {1} and number of '
                             'unpaired electrons = {2}'.format(self.to_adjacency_list(), multiplicity, n_rad))
        if raise_charge_exception:
            if self.get_net_charge() != 0:
                raise ValueError('Non-neutral molecule encountered. '
                                 'Currently, RMG does not support ion chemistry.\n {0}'.format(adjlist))
        return self

    def from_xyz(self, atomic_nums, coordinates, critical_distance_factor=0.45, raise_atomtype_exception=True):
        """
        Create an RMG molecule from a list of coordinates and a corresponding
        list of atomic numbers. These are typically received from CCLib and the
        molecule is sent to `ConnectTheDots` so will only contain single bonds.
        """

        _rdkit_periodic_table = elements.GetPeriodicTable()

        for i, at_num in enumerate(atomic_nums):
            atom = Atom(_rdkit_periodic_table.GetElementSymbol(int(at_num)))
            atom.coords = coordinates[i]
            self.add_atom(atom)
        return self.connect_the_dots(critical_distance_factor=critical_distance_factor, raise_atomtype_exception=raise_atomtype_exception)

    def to_single_bonds(self, raise_atomtype_exception=True):
        """
        Returns a copy of the current molecule, consisting of only single bonds.
        
        This is useful for isomorphism comparison against something that was made
        via from_xyz, which does not attempt to perceive bond orders
        """
        cython.declare(atom1=Atom, atom2=Atom, bond=Bond, newMol=Molecule, atoms=list, mapping=dict)

        new_mol = Molecule()
        atoms = self.atoms
        mapping = {}
        for atom1 in atoms:
            atom2 = new_mol.add_atom(Atom(atom1.element))
            mapping[atom1] = atom2

        for atom1 in atoms:
            for atom2 in atom1.bonds:
                bond = Bond(mapping[atom1], mapping[atom2], 1)
                new_mol.add_bond(bond)
        new_mol.update_atomtypes(raise_exception=raise_atomtype_exception)
        return new_mol

    def to_inchi(self, backend='rdkit-first'):
        """
        Convert a molecular structure to an InChI string. Uses
        `RDKit <https://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity.

        or

        Convert a molecular structure to an InChI string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.

        It is possible to use a single backend or try different backends in sequence.
        The available options for the ``backend`` argument: 'rdkit-first'(default),
        'openbabel-first', 'rdkit', or 'openbabel'.
        """
        try:
            return translator.to_inchi(self, backend=backend)
        except:
            logging.exception(f"Error for molecule \n{self.to_adjacency_list()}")
            raise

    def to_augmented_inchi(self, backend='rdkit-first'):
        """
        Adds an extra layer to the InChI denoting the multiplicity
        of the molecule.

        Separate layer with a forward slash character.

        RDKit and Open Babel are the two backends used in RMG. It is possible to use a
        single backend or try different backends in sequence. The available options for the ``backend``
        argument: 'rdkit-first'(default), 'openbabel-first', 'rdkit', or 'openbabel'.
        """
        try:
            return translator.to_inchi(self, backend=backend, aug_level=2)
        except:
            logging.exception(f"Error for molecule \n{self.to_adjacency_list()}")
            raise

    def to_inchi_key(self, backend='rdkit-first'):
        """
        Convert a molecular structure to an InChI Key string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.

        or

        Convert a molecular structure to an InChI Key string. Uses
        `RDKit <https://rdkit.org/>`_ to perform the conversion.

        It is possible to use a single backend or try different backends in sequence.
        The available options for the ``backend`` argument: 'rdkit-first'(default),
        'openbabel-first', 'rdkit', or 'openbabel'.
        """
        try:
            return translator.to_inchi_key(self, backend=backend)
        except:
            logging.exception(f"Error for molecule \n{self.to_adjacency_list()}")
            raise

    def to_augmented_inchi_key(self, backend='rdkit-first'):
        """
        Adds an extra layer to the InChIKey denoting the multiplicity
        of the molecule.

        Simply append the multiplicity string, do not separate by a
        character like forward slash.

        RDKit and Open Babel are the two backends used in RMG. It is possible to use a
        single backend or try different backends in sequence. The available options for the ``backend``
        argument: 'rdkit-first'(default), 'openbabel-first', 'rdkit', or 'openbabel'.
        """
        try:
            return translator.to_inchi_key(self, backend=backend, aug_level=2)
        except:
            logging.exception(f"Error for molecule \n{self.to_adjacency_list()}")
            raise

    def to_smarts(self):
        """
        Convert a molecular structure to an SMARTS string. Uses
        `RDKit <https://rdkit.org/>`_ to perform the conversion.
        Perceives aromaticity and removes Hydrogen atoms.
        """
        return translator.to_smarts(self)

    def to_smiles(self):
        """
        Convert a molecular structure to an SMILES string. 
        
        If there is a Nitrogen atom present it uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion,
        and the SMILES may or may not be canonical.
        
        Otherwise, it uses `RDKit <https://rdkit.org/>`_ to perform the 
        conversion, so it will be canonical SMILES.
        While converting to an RDMolecule it will perceive aromaticity
        and removes Hydrogen atoms.
        """

        return translator.to_smiles(self)

    def to_rdkit_mol(self, *args, **kwargs):
        """
        Convert a molecular structure to a RDKit rdmol object.
        """
        return converter.to_rdkit_mol(self, *args, **kwargs)

    def to_adjacency_list(self, label='', remove_h=False, remove_lone_pairs=False, old_style=False):
        """
        Convert the molecular structure to a string adjacency list.
        """
        from rmgpy.molecule.adjlist import to_adjacency_list
        result = to_adjacency_list(self.vertices, self.multiplicity, metal=self.metal, facet=self.facet, 
                                   label=label, group=False, remove_h=remove_h,
                                   remove_lone_pairs=remove_lone_pairs, old_style=old_style)
        return result

    def find_h_bonds(self):
        """
        generates a list of (new-existing H bonds ignored) possible Hbond coordinates [(i1,j1),(i2,j2),...] where i and j values
        correspond to the indexes of the atoms involved, Hbonds are allowed if they meet
        the following constraints:

           1) between a H and [O,N] atoms
           2) the hydrogen is covalently bonded to an O or N
           3) the Hydrogen bond must complete a ring with at least 5 members
           4) An atom can only be hydrogen bonded to one other atom
        """
        pot_bonds = []

        ONatoms = [a for a in self.atoms if a.is_oxygen() or a.is_nitrogen()]
        ONinds = [n for n, a in enumerate(self.atoms) if a.is_oxygen() or a.is_nitrogen()]

        for i, atm1 in enumerate(self.atoms):
            if atm1.atomtype.label == 'H0':
                atm_covs = [q for q in atm1.bonds.keys()]
                if len(atm_covs) > 1:  # H is already H bonded
                    continue
                else:
                    atm_cov = atm_covs[0]
                if (atm_cov.is_oxygen() or atm_cov.is_nitrogen()):  # this H can be H-bonded
                    for k, atm2 in enumerate(ONatoms):
                        if all([not np.isclose(0.1, q.order) for q in
                                atm2.bonds.values()]):  # atm2 not already H bonded
                            dist = len(find_shortest_path(atm1, atm2)) - 1
                            if dist > 3:
                                j = ONinds[k]
                                pot_bonds.append((i, j))
        return pot_bonds

    def generate_h_bonded_structures(self):
        """
        generates a list of Hbonded molecular structures in addition to the
        constraints on Hydrogen bonds applied in the find_H_Bonds function
        the generated structures are constrained to:

            1) An atom can only be hydrogen bonded to one other atom
            2) Only two H-bonds can exist in a given molecule

        the second is done to avoid explosive growth in the number of 
        structures as without this constraint the number of possible 
        structures grows 2^n where n is the number of possible H-bonds
        """
        structs = []
        Hbonds = self.find_h_bonds()
        for i, bd1 in enumerate(Hbonds):
            molc = self.copy(deep=True)
            molc.add_bond(Bond(molc.atoms[bd1[0]], molc.atoms[bd1[1]], order=0.1))
            structs.append(molc)
            for j, bd2 in enumerate(Hbonds):
                if j < i and bd1[0] != bd2[0] and bd1[1] != bd2[1]:
                    molc = self.copy(deep=True)
                    molc.add_bond(Bond(molc.atoms[bd1[0]], molc.atoms[bd1[1]], order=0.1))
                    molc.add_bond(Bond(molc.atoms[bd2[0]], molc.atoms[bd2[1]], order=0.1))
                    structs.append(molc)

        return structs

    def remove_h_bonds(self):
        """
        removes any present hydrogen bonds from the molecule
        """

        atoms = self.atoms
        for i, atm1 in enumerate(atoms):
            for j, atm2 in enumerate(atoms):
                if j < i and self.has_bond(atm1, atm2):
                    bd = self.get_bond(atm1, atm2)
                    if np.isclose(0.1, bd.order):
                        self.remove_bond(bd)
        return

    def is_linear(self):
        """
        Return :data:`True` if the structure is linear and :data:`False`
        otherwise.
        """

        atom_count = len(self.vertices)

        # Monatomic molecules are definitely nonlinear
        if atom_count == 1:
            return False
        # Diatomic molecules are definitely linear
        elif atom_count == 2:
            return True
        # Cyclic molecules are definitely nonlinear
        elif self.is_cyclic():
            return False

        # True if all bonds are double bonds (e.g. O=C=O)
        all_double_bonds = True
        for atom1 in self.vertices:
            for bond in atom1.edges.values():
                if not bond.is_double(): all_double_bonds = False
        if all_double_bonds: return True

        # True if alternating single-triple bonds (e.g. H-C#C-H)
        # This test requires explicit hydrogen atoms
        for atom in self.vertices:
            bonds = list(atom.edges.values())
            if len(bonds) == 1:
                continue  # ok, next atom
            if len(bonds) > 2:
                break  # fail!
            if bonds[0].is_single() and bonds[1].is_triple():
                continue  # ok, next atom
            if bonds[1].is_single() and bonds[0].is_triple():
                continue  # ok, next atom
            break  # fail if we haven't continued
        else:
            # didn't fail
            return True

        # not returned yet? must be nonlinear
        return False

    def is_aromatic(self):
        """ 
        Returns ``True`` if the molecule is aromatic, or ``False`` if not.  
        Iterates over the SSSR's and searches for rings that consist solely of Cb 
        atoms.  Assumes that aromatic rings always consist of 6 atoms. 
        In cases of naphthalene, where a 6 + 4 aromatic system exists,
        there will be at least one 6 membered aromatic ring so this algorithm
        will not fail for fused aromatic rings.
        """
        cython.declare(rc=list, cycle=list, atom=Atom)
        rc = self.get_relevant_cycles()
        if rc:
            for cycle in rc:
                if len(cycle) == 6:
                    for atom in cycle:
                        # print atom.atomtype.label
                        if atom.atomtype.label == 'Cb' or atom.atomtype.label == 'Cbf':
                            continue
                            # Go onto next cycle if a non Cb atomtype was discovered in this cycle
                        break
                    else:
                        # Molecule is aromatic when all 6 atoms are type 'Cb'
                        return True
        return False

    def is_heterocyclic(self):
        """
        Returns ``True`` if the molecule is heterocyclic, or ``False`` if not.
        """
        if self.is_cyclic():
            for atom in self.atoms:
                if atom.is_non_hydrogen() and not atom.is_carbon() and self.is_vertex_in_cycle(atom):
                    return True
        return False

    def count_internal_rotors(self):
        """
        Determine the number of internal rotors in the structure. Any single
        bond not in a cycle and between two atoms that also have other bonds
        is considered to be a pivot of an internal rotor.
        """
        count = 0
        for atom1 in self.vertices:
            for atom2, bond in atom1.edges.items():
                if (self.vertices.index(atom1) < self.vertices.index(atom2) and
                        bond.is_single() and not self.is_bond_in_cycle(bond)):
                    if len(atom1.edges) > 1 and len(atom2.edges) > 1:
                        count += 1
        return count

    def calculate_cp0(self):
        """
        Return the value of the heat capacity at zero temperature in J/mol*K.
        """
        if self.contains_surface_site():
            return 0.01
        if len(self.atoms) == 1:
            return 2.5 * constants.R
        else:
            return (3.5 if self.is_linear() else 4.0) * constants.R

    def calculate_cpinf(self):
        """
        Return the value of the heat capacity at infinite temperature in J/mol*K.
        """
        cython.declare(n_atoms=cython.int, n_vib=cython.int, n_rotors=cython.int)

        if self.contains_surface_site():
            # ToDo: internal rotors could still act as rotors
            return constants.R * 3 * len(self.vertices)

        if len(self.vertices) == 1:
            return self.calculate_cp0()
        else:
            n_atoms = len(self.vertices)
            n_vib = 3 * n_atoms - (5 if self.is_linear() else 6)
            n_rotors = self.count_internal_rotors()
            n_vib -= n_rotors

            return self.calculate_cp0() + (n_vib + 0.5 * n_rotors) * constants.R

    def get_symmetry_number(self):
        """
        Returns the symmetry number of Molecule.
        First checks whether the value is stored as an attribute of Molecule.
        If not, it calls the calculate_symmetry_number method.
        """
        if self.symmetry_number == -1:
            self.calculate_symmetry_number()
        return self.symmetry_number

    def calculate_symmetry_number(self):
        """
        Return the symmetry number for the structure. The symmetry number
        includes both external and internal modes.
        """
        from rmgpy.molecule.symmetry import calculate_symmetry_number
        self.update_connectivity_values()  # for consistent results
        self.symmetry_number = calculate_symmetry_number(self)
        return self.symmetry_number

    def is_radical(self):
        """
        Return ``True`` if the molecule contains at least one radical electron,
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.radical_electrons > 0:
                return True
        return False

    def has_lone_pairs(self):
        """
        Return ``True`` if the molecule contains at least one lone electron pair,
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.lone_pairs > 0:
                return True
        return False
    
    def has_charge(self):
        for atom in self.vertices:
            if atom.charge != 0:
                return True
        return False

    def has_halogen(self):
        """
        Return ``True`` if the molecule contains at least one halogen (F, Cl, Br, or I),
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        for atom in self.vertices:
            if atom.is_halogen():
                return True
        return False

    def is_aryl_radical(self, aromatic_rings=None, save_order=False):
        """
        Return ``True`` if the molecule only contains aryl radicals,
        i.e., radical on an aromatic ring, or ``False`` otherwise.
        If no ``aromatic_rings`` provided, aromatic rings will be searched in-place,
        and this process may involve atom order change by default. Set ``save_order`` to
        ``True`` to force the atom order unchanged.
        """
        cython.declare(atom=Atom, total=int, aromatic_atoms=set, aryl=int)
        if aromatic_rings is None:
            aromatic_rings = self.get_aromatic_rings(save_order=save_order)[0]

        total = self.get_radical_count()
        aromatic_atoms = set([atom for atom in itertools.chain.from_iterable(aromatic_rings)])
        aryl = sum([atom.radical_electrons for atom in aromatic_atoms])

        return total == aryl

    def generate_resonance_structures(self, keep_isomorphic=False, filter_structures=True, save_order=False):
        """Returns a list of resonance structures of the molecule."""

        try:
            return resonance.generate_resonance_structures(self, keep_isomorphic=keep_isomorphic,
                                                       filter_structures=filter_structures,
                                                       save_order=save_order,
                                                       )
        except:
            logging.warning("Resonance structure generation failed for {}".format(self))
            return [self.copy(deep=True)]

    def get_url(self):
        """
        Get a URL to the molecule's info page on the RMG website.
        """
        # eg. http://dev.rmg.mit.edu/database/kinetics/reaction/reactant1=1%20C%200%20%7B2,S%7D;2%20O%200%20%7B1,S%7D;__reactant2=1%20C%202T;__product1=1%20C%201;__product2=1%20C%200%20%7B2,S%7D;2%20O%201%20%7B1,S%7D;

        base_url = "https://rmg.mit.edu/database/molecule/"
        adjlist = self.to_adjacency_list(remove_h=False)
        url = base_url + quote(adjlist)
        return url.strip('_')

    def get_radical_atoms(self):
        """
        Return the atoms in the molecule that have unpaired electrons.
        """
        radical_atoms_list = []
        for atom in self.vertices:
            if atom.radical_electrons > 0:
                radical_atoms_list.append(atom)
        return radical_atoms_list

    def update_lone_pairs(self):
        """
        Iterate through the atoms in the structure and calculate the
        number of lone electron pairs, assuming a neutral molecule.
        """
        cython.declare(atom1=Atom, atom2=Atom, bond12=Bond, order=float)
        for atom1 in self.vertices:
            if atom1.is_hydrogen() or atom1.is_surface_site() or atom1.is_electron() or atom1.is_lithium():
                atom1.lone_pairs = 0
            else:
                order = atom1.get_total_bond_order()
                atom1.lone_pairs = (elements.PeriodicSystem.valence_electrons[atom1.symbol]
                                   - atom1.radical_electrons - atom1.charge - int(order)) / 2.0
                if atom1.lone_pairs % 1 > 0 or atom1.lone_pairs > 4:
                    logging.error("Unable to determine the number of lone pairs for "
                                  "element {0} in {1}".format(atom1, self))

    def get_net_charge(self):
        """
        Iterate through the atoms in the structure and calculate the net charge
        on the overall molecule.
        """
        return sum([atom.charge for atom in self.vertices])

    def get_charge_span(self):
        """
        Iterate through the atoms in the structure and calculate the charge span
        on the overall molecule.
        The charge span is a measure of the number of charge separations in a molecule.
        """
        abs_net_charge = abs(self.get_net_charge())
        sum_of_abs_charges = sum([abs(atom.charge) for atom in self.vertices])
        return (sum_of_abs_charges - abs_net_charge) / 2

    def saturate_unfilled_valence(self, update=True):
        """
        Saturate the molecule by adding H atoms to any unfilled valence
        """

        saturator = Saturator()
        saturator.saturate(self.atoms)
        if update: self.update()

    def saturate_radicals(self, raise_atomtype_exception=True):
        """
        Saturate the molecule by replacing all radicals with bonds to hydrogen atoms.  Changes self molecule object.
        """
        cython.declare(added=dict, atom=Atom, i=int, H=Atom, bond=Bond)
        added = {}
        for atom in self.atoms:
            for i in range(atom.radical_electrons):
                H = Atom('H', radical_electrons=0, lone_pairs=0, charge=0)
                bond = Bond(atom, H, 1)
                self.add_atom(H)
                self.add_bond(bond)
                if atom not in added:
                    added[atom] = []
                added[atom].append([H, bond])
                atom.decrement_radical()

        # Update the atom types of the saturated structure (not sure why
        # this is necessary, because saturating with H shouldn't be
        # changing atom types, but it doesn't hurt anything and is not
        # very expensive, so will do it anyway)
        self.sort_atoms()
        self.update_atomtypes(raise_exception=raise_atomtype_exception)
        self.multiplicity = 1

        return added

    def replace_halogen_with_hydrogen(self, raise_atomtype_exception=True):
        """
        Replace all halogens in a molecule with hydrogen atoms. Changes self molecule object.
        """
        cython.declare(halogen_atom_list=list, atom=Atom, bond_to_replace_dict=dict, H_atom=Atom,
                       bonded_atom=Atom, bond_to_replace=Bond, new_bond=Bond)
        # the list of halogen atoms must be obtained before any of the halogen atoms are replaced because it changes
        # the order of self.atoms
        halogen_atom_list = [atom for atom in self.atoms if atom.is_halogen()]
        for atom in halogen_atom_list:
            if not atom.charge == 0:
                raise ValueError('For a given molecule {0}, a halogen atom {1} with charge {2} cannot be replaced '
                                 'with a hydrogen atom'.format(self.to_smiles(), atom.symbol, atom.charge))
            bond_to_replace_dict = self.get_bonds(atom)
            self.remove_atom(atom)
            H_atom = Atom('H', radical_electrons=atom.radical_electrons, lone_pairs=0, charge=0)
            self.add_atom(H_atom)
            for bonded_atom, bond_to_replace in bond_to_replace_dict.items():
                new_bond = Bond(H_atom, bonded_atom, order=bond_to_replace.order)
                self.add_bond(new_bond)

        # Update the atom types of the new structure
        self.sort_atoms()
        self.update_atomtypes(raise_exception=raise_atomtype_exception)

    def to_group(self):
        """
        This method converts a list of atoms in a Molecule to a Group object.
        """

        # Create GroupAtom object for each atom in the molecule
        group_atoms = OrderedDict()  # preserver order of atoms in original container
        for atom in self.atoms:
            group_atoms[atom] = gr.GroupAtom(atomtype=[atom.atomtype],
                                             radical_electrons=[atom.radical_electrons],
                                             charge=[atom.charge],
                                             lone_pairs=[atom.lone_pairs],
                                             label=atom.label,
                                             site=[atom.site] if atom.site else [],
                                             morphology=[atom.morphology] if atom.morphology else [],
                                             )

        group = gr.Group(atoms=list(group_atoms.values()), multiplicity=[self.multiplicity], metal=[self.metal] if self.metal else [],
                         facet=[self.facet] if self.facet else [])

        # Create GroupBond for each bond between atoms in the molecule
        for atom in self.atoms:
            for bonded_atom, bond in atom.edges.items():
                group.add_bond(gr.GroupBond(group_atoms[atom], group_atoms[bonded_atom], order=[bond.order]))

        group.update()

        return group

    def identify_ring_membership(self):
        """
        Performs ring perception and saves ring membership information to the Atom.props attribute.
        """
        cython.declare(rc=list, atom=Atom, ring=list)

        # Get the set of relevant cycles
        rc = self.get_relevant_cycles()
        # Identify whether each atom is in a ring
        for atom in self.atoms:
            atom.props['inRing'] = False
            for ring in rc:
                if atom in ring:
                    atom.props['inRing'] = True
                    break

    def count_aromatic_rings(self):
        """
        Count the number of aromatic rings in the current molecule, as determined by the benzene bond type.
        This is purely dependent on representation and is unrelated to the actual aromaticity of the molecule.

        Returns an integer corresponding to the number or aromatic rings.
        """
        cython.declare(rings=list, count=int, ring=list, bonds=list, bond=Bond)
        rings = self.get_relevant_cycles()
        count = 0
        for ring in rings:
            if len(ring) != 6:
                # We currently only recognize 6-membered rings as being aromatic
                continue

            bonds = self.get_edges_in_cycle(ring)
            if all([bond.is_benzene() for bond in bonds]):
                count += 1

        return count

    def get_aromatic_rings(self, rings=None, save_order=False):
        """
        Returns all aromatic rings as a list of atoms and a list of bonds.

        Identifies rings using `Graph.get_smallest_set_of_smallest_rings()`, then uses RDKit to perceive aromaticity.
        RDKit uses an atom-based pi-electron counting algorithm to check aromaticity based on Huckel's Rule.
        Therefore, this method identifies "true" aromaticity, rather than simply the RMG bond type.

        The method currently restricts aromaticity to six-membered carbon-only rings. This is a limitation imposed
        by RMG, and not by RDKit.

        By default, the atom order will be sorted to get consistent results from different runs. The atom order can
        be saved when dealing with problems that are sensitive to the atom map.
        """
        cython.declare(rd_atom_indices=dict, ob_atom_ids=dict, aromatic_rings=list, aromatic_bonds=list)
        cython.declare(ring0=list, i=cython.int, atom1=Atom, atom2=Atom)

        from rdkit.Chem.rdchem import BondType
        AROMATIC = BondType.AROMATIC

        if rings is None:
            rings = self.get_relevant_cycles()

        # Remove rings that share more than 3 atoms, since they cannot be planar
        cython.declare(toRemove=set, j=cython.int, toRemoveSorted=list)
        if len(rings) < 2:
            pass
        else:
            to_remove = set()
            for i, j in itertools.combinations(range(len(rings)), 2):
                if len(set(rings[i]) & set(rings[j])) > 2:
                    to_remove.add(i)
                    to_remove.add(j)

            to_remove_sorted = sorted(to_remove, reverse=True)

            for i in to_remove_sorted:
                del rings[i]

        # Only keep rings with exactly 6 atoms, since RMG can only handle aromatic benzene
        rings = [ring for ring in rings if len(ring) == 6]

        if not rings:
            return [], []

        try:
            rdkitmol, rd_atom_indices = converter.to_rdkit_mol(self, remove_h=False, return_mapping=True, save_order=save_order)
        except ValueError:
            logging.warning('Unable to check aromaticity by converting to RDKit Mol.')
        else:
            aromatic_rings = []
            aromatic_bonds = []
            for ring0 in rings:
                aromatic_bonds_in_ring = []
                # Figure out which atoms and bonds are aromatic and reassign appropriately:
                for i, atom1 in enumerate(ring0):
                    if not atom1.is_carbon():
                        # all atoms in the ring must be carbon in RMG for our definition of aromatic
                        break
                    for atom2 in ring0[i + 1:]:
                        if self.has_bond(atom1, atom2):
                            # Check for aromaticity using the bond type rather than GetIsAromatic because
                            # aryne triple bonds return True for GetIsAromatic but are not aromatic bonds
                            if rdkitmol.GetBondBetweenAtoms(rd_atom_indices[atom1],
                                                            rd_atom_indices[atom2]).GetBondType() is AROMATIC:
                                aromatic_bonds_in_ring.append(self.get_bond(atom1, atom2))
                else:  # didn't break so all atoms are carbon
                    if len(aromatic_bonds_in_ring) == 6:
                        aromatic_rings.append(ring0)
                        aromatic_bonds.append(aromatic_bonds_in_ring)

            return aromatic_rings, aromatic_bonds

        logging.info('Trying to use OpenBabel to check aromaticity.')
        try:
            obmol, ob_atom_ids = converter.to_ob_mol(self, return_mapping=True, save_order=save_order)
        except DependencyError:
            logging.warning('Unable to check aromaticity by converting for OB Mol.')
            return [], []
        else:
            aromatic_rings = []
            aromatic_bonds = []
            for ring0 in rings:
                aromatic_bonds_in_ring = []
                # Figure out which atoms and bonds are aromatic and reassign appropriately:
                for i, atom1 in enumerate(ring0):
                    if not atom1.is_carbon():
                        # all atoms in the ring must be carbon in RMG for our definition of aromatic
                        break
                    for atom2 in ring0[i + 1:]:
                        if self.has_bond(atom1, atom2):
                            if obmol.GetBond(obmol.GetAtomById(ob_atom_ids[atom1]),
                                             obmol.GetAtomById(ob_atom_ids[atom2])).IsAromatic():
                                aromatic_bonds_in_ring.append(self.get_bond(atom1, atom2))
                else:  # didn't break so all atoms are carbon
                    if len(aromatic_bonds_in_ring) == 6:
                        aromatic_rings.append(ring0)
                        aromatic_bonds.append(aromatic_bonds_in_ring)

            return aromatic_rings, aromatic_bonds

    def get_deterministic_sssr(self):
        """
        Modified `Graph` method `get_smallest_set_of_smallest_rings` by sorting calculated cycles
        by short length and then high atomic number instead of just short length (for cases where
        multiple cycles with same length are found, `get_smallest_set_of_smallest_rings` outputs
        non-determinstically).
        
        For instance, molecule with this smiles: C1CC2C3CSC(CO3)C2C1, will have non-deterministic
        output from `get_smallest_set_of_smallest_rings`, which leads to non-deterministic bicyclic decomposition.
        Using this new method can effectively prevent this situation.

        Important Note: This method returns an incorrect set of SSSR in certain molecules (such as cubane).
        It is recommended to use the main `Graph.get_smallest_set_of_smallest_rings` method in new applications.
        Alternatively, consider using `Graph.get_relevant_cycles` for deterministic output.

        In future development, this method should ideally be replaced by some method to select a deterministic
        set of SSSR from the set of Relevant Cycles, as that would be a more robust solution.
        """
        cython.declare(vertices=list, vertices_to_remove=list, root_candidates_tups=list, graphs=list)
        cython.declare(cycle_list=list, cycle_candidate_tups=list, cycles=list, cycle0=list, origin_conn_dict=dict)

        cython.declare(graph=Molecule, graph0=Molecule, vertex=Atom, root_vertex=Atom)

        # Make a copy of the graph so we don't modify the original
        graph = self.copy(deep=True)
        vertices = graph.vertices[:]

        # Step 1: Remove all terminal vertices
        done = False
        while not done:
            vertices_to_remove = []
            for vertex in graph.vertices:
                if len(vertex.edges) == 1: vertices_to_remove.append(vertex)
            done = len(vertices_to_remove) == 0
            # Remove identified vertices from graph
            for vertex in vertices_to_remove:
                graph.remove_vertex(vertex)

        graph.update_connectivity_values()
        # get original connectivity values
        origin_conn_dict = {}
        for v in graph.vertices:
            origin_conn_dict[v] = get_vertex_connectivity_value(v)

        # Step 2: Remove all other vertices that are not part of cycles
        vertices_to_remove = []
        for vertex in graph.vertices:
            found = graph.is_vertex_in_cycle(vertex)
            if not found:
                vertices_to_remove.append(vertex)
        # Remove identified vertices from graph
        for vertex in vertices_to_remove:
            graph.remove_vertex(vertex)

        # Step 3: Split graph into remaining subgraphs
        graphs = graph.split()

        # Step 4: Find ring sets in each subgraph
        cycle_list = []
        for graph0 in graphs:

            while len(graph0.vertices) > 0:

                # Choose root vertex as vertex with smallest number of edges
                root_vertex = None
                graph0.update_connectivity_values()

                root_candidates_tups = []
                for vertex in graph0.vertices:
                    tup = (vertex, get_vertex_connectivity_value(vertex), -origin_conn_dict[vertex])
                    root_candidates_tups.append(tup)

                root_vertex = sorted(root_candidates_tups, key=_skip_first, reverse=True)[0][0]

                # Get all cycles involving the root vertex
                cycles = graph0.get_all_cycles(root_vertex)
                if len(cycles) == 0:
                    # This vertex is no longer in a ring, so remove it
                    graph0.remove_vertex(root_vertex)
                    continue

                # Keep the smallest of the cycles found above
                cycle_candidate_tups = []
                for cycle0 in cycles:
                    tup = (cycle0, len(cycle0), -sum([origin_conn_dict[v] for v in cycle0]),
                           -sum([v.element.number for v in cycle0]),
                           -sum([v.get_total_bond_order() for v in cycle0]))
                    cycle_candidate_tups.append(tup)

                cycle = sorted(cycle_candidate_tups, key=_skip_first)[0][0]

                cycle_list.append(cycle)

                # Remove the root vertex to create single edges, note this will not
                # function properly if there is no vertex with 2 edges (i.e. cubane)
                graph0.remove_vertex(root_vertex)

                # Remove from the graph all vertices in the cycle that have only one edge
                lone_carbon = True
                while lone_carbon:
                    lone_carbon = False
                    vertices_to_remove = []

                    for vertex in cycle:
                        if len(vertex.edges) == 1:
                            lone_carbon = True
                            vertices_to_remove.append(vertex)
                    else:
                        for vertex in vertices_to_remove:
                            graph0.remove_vertex(vertex)

        # Map atoms in cycles back to atoms in original graph
        for i in range(len(cycle_list)):
            cycle_list[i] = [self.vertices[vertices.index(v)] for v in cycle_list[i]]

        return cycle_list

    def kekulize(self):
        """
        Kekulizes an aromatic molecule.
        """
        kekulize(self)

    def assign_atom_ids(self):
        """
        Assigns an index to every atom in the molecule for tracking purposes.
        Uses entire range of cython's integer values to reduce chance of duplicates
        """

        global atom_id_counter

        for atom in self.atoms:
            atom.id = atom_id_counter
            atom_id_counter += 1
            if atom_id_counter == 2 ** 15:
                atom_id_counter = -2 ** 15

    def atom_ids_valid(self):
        """
        Checks to see if the atom IDs are valid in this structure
        """
        num_atoms = len(self.atoms)
        num_ids = len(set([atom.id for atom in self.atoms]))

        if num_atoms == num_ids:
            # all are unique
            return True
        return False

    def is_identical(self, other, strict=True):
        """
        Performs isomorphism checking, with the added constraint that atom IDs must match.

        Primary use case is tracking atoms in reactions for reaction degeneracy determination.

        Returns :data:`True` if two graphs are identical and :data:`False` otherwise.

        If ``strict=False``, performs the check ignoring electrons and resonance structures.
        """
        cython.declare(atom_ids=set, other_ids=set, atom_list=list, other_list=list, mapping=dict)
        from rmgpy.molecule.fragment import Fragment

        if not isinstance(other, (Molecule, Fragment)):
            raise TypeError(
                'Got a {0} object for parameter "other", when a Molecule object is required.'.format(other.__class__))

        # Get a set of atom indices for each molecule
        atom_ids = set([atom.id for atom in self.atoms])
        other_ids = set([atom.id for atom in other.atoms])

        if atom_ids == other_ids:
            # If the two molecules have the same indices, then they might be identical
            # Sort the atoms by ID
            atom_list = sorted(self.atoms, key=attrgetter('id'))
            other_list = sorted(other.atoms, key=attrgetter('id'))

            # If matching atom indices gives a valid mapping, then the molecules are fully identical
            mapping = {}
            for atom1, atom2 in zip(atom_list, other_list):
                mapping[atom1] = atom2

            return self.is_mapping_valid(other, mapping, equivalent=True, strict=strict)
        else:
            # The molecules don't have the same set of indices, so they are not identical
            return False

    def get_nth_neighbor(self, starting_atoms, distance_list, ignore_list=None, n=1):
        """
        Recursively get the Nth nonHydrogen neighbors of the starting_atoms, and return them in a list.
        `starting_atoms` is a list of :class:Atom for which we will get the nth neighbor.
        `distance_list` is a list of integers, corresponding to the desired neighbor distances.
        `ignore_list` is a list of :class:Atom that have been counted in (n-1)th neighbor, and will not be returned.
        `n` is an integer, corresponding to the distance to be calculated in the current iteration.
        """
        if ignore_list is None:
            ignore_list = []

        neighbors = []
        for atom in starting_atoms:
            new_neighbors = [neighbor for neighbor in self.get_bonds(atom) if neighbor.is_non_hydrogen()]
            neighbors.extend(new_neighbors)

        neighbors = list(set(neighbors) - set(ignore_list))
        for atom in starting_atoms:
            ignore_list.append(atom)
        if n < max(distance_list):
            if n in distance_list:
                neighbors += self.get_nth_neighbor(neighbors, distance_list, ignore_list, n + 1)
            else:
                neighbors = self.get_nth_neighbor(neighbors, distance_list, ignore_list, n + 1)
        return neighbors

    def enumerate_bonds(self):
        """
        Count the number of each type of bond (e.g. 'C-H', 'C=C') present in the molecule
        :return: dictionary, with bond strings as keys and counts as values
        """
        bond_count = defaultdict(int)
        bonds = self.get_all_edges()

        for bond in bonds:
            bond_count[bond.get_bond_string()] += 1

        return dict(bond_count)

    def get_surface_sites(self):
        """
        Get a list of surface site atoms in the molecule.
        Returns:
            List(Atom): A list containing the surface site atoms in the molecule
        """
        cython.declare(atom=Atom)
        return [atom for atom in self.atoms if atom.is_surface_site()]

    def is_multidentate(self):
        """
        Return ``True`` if the adsorbate contains at least two binding sites,
        or ``False`` otherwise.
        """
        cython.declare(atom=Atom)
        if len([atom for atom in self.atoms if atom.is_surface_site()])>=2:
            return True
        return False

    def get_adatoms(self):
        """
        Get a list of adatoms in the molecule.
        Returns:
            List(Atom): A list containing the adatoms in the molecule
        """
        cython.declare(surface_site=Atom, adatoms=list)
        adatoms = []
        for surface_site in self.get_surface_sites():
            if surface_site.bonds:
                adatoms.extend(surface_site.bonds.keys())
        return adatoms

    def get_desorbed_molecules(self):
        """
        Get a list of desorbed molecules by desorbing the molecule from the surface.
        
        Returns a list of Molecules.  Each molecule's atoms will be labeled corresponding to
        the bond order with the surface:
        ``*1`` - single bond
        ``*2`` - double bond
        ``*3`` - triple bond
        ``*4`` - quadruple bond
        """
        cython.declare(desorbed_molecules=list, desorbed_molecule=Molecule, sites_to_remove=list, adsorbed_atoms=list,
                       site=Atom, numbonds=cython.int, bonded_atom=Atom, bond=Bond, i=cython.int, j=cython.int, atom0=Atom,
                       atom1=Atom)

        if not self.contains_surface_site():
            return []

        desorbed_molecule = self.copy(deep=True)
        desorbed_molecule.clear_labeled_atoms()
        sites_to_remove = desorbed_molecule.get_surface_sites()
        adsorbed_atoms = []
        for site in sites_to_remove:
            numbonds = len(site.bonds)
            if numbonds == 0:
                # vanDerWaals
                pass
            else:
                assert len(site.bonds) == 1, "Each surface site can only be bonded to 1 atom"
                (bonded_atom, bond), = site.bonds.items()
                adsorbed_atoms.append(bonded_atom)
                desorbed_molecule.remove_bond(bond)
                if bond.is_single():
                    bonded_atom.increment_radical()
                    bonded_atom.label = '*1'
                elif bond.is_double():
                    bonded_atom.increment_radical()
                    bonded_atom.increment_radical()
                    bonded_atom.label = '*2'
                elif bond.is_triple():
                    bonded_atom.increment_radical()
                    bonded_atom.increment_lone_pairs()
                    bonded_atom.label = '*3'
                elif bond.is_quadruple():
                    bonded_atom.increment_radical()
                    bonded_atom.increment_radical()
                    bonded_atom.increment_lone_pairs()
                    bonded_atom.label = '*4'
                else:
                    raise NotImplementedError("Can't remove surface bond of type {}".format(bond.order))
            desorbed_molecule.remove_atom(site)

        desorbed_molecules = [desorbed_molecule.copy(deep=True)]
        if len(adsorbed_atoms) > 1:
            # multidentate adsorption.
            # Try to turn adjacent biradical into a bond.
            for i,j in itertools.combinations(range(len(adsorbed_atoms)),2):
                try:
                    atom0 = adsorbed_atoms[i]
                    atom1 = adsorbed_atoms[j]
                    bond = atom0.bonds[atom1]
                except KeyError:
                    pass # the two adsorbed atoms are not bonded to each other
                else:
                    if (atom0.radical_electrons and atom1.radical_electrons and bond.order < 3):
                        bond.increment_order()
                        atom0.decrement_radical()
                        atom1.decrement_radical()
                        desorbed_molecules.append(desorbed_molecule.copy(deep=True))
                        if (atom0.radical_electrons and
                                atom1.radical_electrons and
                                bond.order < 3):
                            # There are still spare adjacenct radicals, so do it again
                            bond.increment_order()
                            atom0.decrement_radical()
                            atom1.decrement_radical()
                            desorbed_molecules.append(desorbed_molecule.copy(deep=True))
                        if (atom0.lone_pairs and
                                atom1.lone_pairs and 
                                bond.order < 3):
                            # X#C-C#X will end up with .:C-C:. in gas phase
                            # and we want to get to .C#C. but not :C=C:
                            bond.increment_order()
                            atom0.decrement_lone_pairs()
                            atom0.increment_radical()
                            atom1.decrement_lone_pairs()
                            atom1.increment_radical()
                            desorbed_molecules.append(desorbed_molecule.copy(deep=True))
                    #For bidentate CO because we want C[-1]#O[+1] but not .C#O.
                    if (bond.order == 3 and atom0.radical_electrons and 
                        atom1.radical_electrons and 
                        (atom0.lone_pairs or atom1.lone_pairs)):
                        atom0.decrement_radical()
                        atom1.decrement_radical()
                        if atom0.lone_pairs:
                            atom1.increment_lone_pairs()
                        else:
                            atom0.increment_lone_pairs()
                        desorbed_molecules.append(desorbed_molecule.copy(deep=True))

        for desorbed_molecule in desorbed_molecules[:]:
            try:
                desorbed_molecule.update_connectivity_values()
                desorbed_molecule.update()
            except AtomTypeError:
                desorbed_molecules.remove(desorbed_molecule)
                logging.debug(f"Removing {desorbed_molecule} from possible structure list:\n{desorbed_molecule.to_adjacency_list()}")
            else:
                logging.debug("After removing from surface:\n" + desorbed_molecule.to_adjacency_list())

        return desorbed_molecules

# this variable is used to name atom IDs so that there are as few conflicts by 
# using the entire space of integer objects
atom_id_counter = -2 ** 15
