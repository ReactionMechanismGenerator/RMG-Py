#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
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
This module provides classes and methods for working with molecular substructure
patterns. These enable molecules to be searched for common motifs (e.g.
reaction sites).

.. _atom-types:

We define the following basic atom types:

    =============== ============================================================
    Atom type       Description
    =============== ============================================================
    *General atom types*
    ----------------------------------------------------------------------------
    ``R``           any atom with any local bond structure
    ``R!H``         any non-hydrogen atom with any local bond structure
    *Carbon atom types*
    ----------------------------------------------------------------------------
    ``C``           carbon atom with any local bond structure
    ``Cs``          carbon atom with four single bonds
    ``Cd``          carbon atom with one double bond (to carbon) and two single bonds
    ``Cdd``         carbon atom with two double bonds
    ``Ct``          carbon atom with one triple bond and one single bond
    ``CO``          carbon atom with one double bond (to oxygen) and two single bonds
    ``Cb``          carbon atom with two benzene bonds and one single bond
    ``Cbf``         carbon atom with three benzene bonds
    *Hydrogen atom types*
    ----------------------------------------------------------------------------
    ``H``           hydrogen atom with one single bond
    *Oxygen atom types*
    ----------------------------------------------------------------------------
    ``O``           oxygen atom with any local bond structure
    ``Os``          oxygen atom with two single bonds
    ``Od``          oxygen atom with one double bond
    ``Oa``          oxygen atom with no bonds
    =============== ============================================================

.. _bond-types:

We define the following bond types:

    =============== ============================================================
    Bond type       Description
    =============== ============================================================
    ``S``           a single bond
    ``D``           a double bond
    ``T``           a triple bond
    ``B``           a benzene bond
    =============== ============================================================

.. _reaction-recipe-actions:

We define the following reaction recipe actions:

    ============= ============================= ================================
    Action name   Arguments                     Action
    ============= ============================= ================================
    CHANGE_BOND   `center1`, `order`, `center2` change the bond order of the bond between `center1` and `center2` by `order`; do not break or form bonds
    FORM_BOND     `center1`, `order`, `center2` form a new bond between `center1` and `center2` of type `order`
    BREAK_BOND    `center1`, `order`, `center2` break the bond between `center1` and `center2`, which should be of type `order`
    GAIN_RADICAL  `center`, `radical`           increase the number of free electrons on `center` by `radical`
    LOSE_RADICAL  `center`, `radical`           decrease the number of free electrons on `center` by `radical`
    ============= ============================= ================================

"""

import cython

import graph
from exception import ChemPyError

################################################################################

def getAtomType(atom, bonds):
    """
    Determine the appropriate atom type for an :class:`Atom` object `atom`
    with local bond structure `bonds`, a ``dict`` containing atom-bond pairs.
    """

    atomType = ''

    # Count numbers of each higher-order bond type
    double = 0; doubleO = 0; triple = 0; benzene = 0
    for atom2, bond12 in bonds.iteritems():
        if bond12.isDouble():
            if atom2.isOxygen(): doubleO +=1
            else:                double += 1
        elif bond12.isTriple(): triple += 1
        elif bond12.isBenzene(): benzene += 1

    # Use element and counts to determine proper atom type
    if atom.symbol == 'C':
        if   double == 0 and doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Cs'
        elif double == 1 and doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Cd'
        elif double + doubleO == 2        and triple == 0 and benzene == 0: atomType = 'Cdd'
        elif double == 0 and doubleO == 0 and triple == 1 and benzene == 0: atomType = 'Ct'
        elif double == 0 and doubleO == 1 and triple == 0 and benzene == 0: atomType = 'CO'
        elif double == 0 and doubleO == 0 and triple == 0 and benzene == 2: atomType = 'Cb'
        elif double == 0 and doubleO == 0 and triple == 0 and benzene == 3: atomType = 'Cbf'
    elif atom.symbol == 'H':
        atomType = 'H'
    elif atom.symbol == 'O':
        if   double + doubleO == 0 and triple == 0 and benzene == 0: atomType = 'Os'
        elif double + doubleO == 1 and triple == 0 and benzene == 0: atomType = 'Od'
        elif len(bonds) == 0:                                        atomType = 'Oa'

    # Raise exception if we could not identify the proper atom type
    if atomType == '':
        raise ChemPyError('Unable to determine atom type for atom %s.' % atom)

    return atomType

def atomTypesEquivalent(atomType1, atomType2):
    """
    Returns ``True`` if two atom types `atomType1` and `atomType2` are
    equivalent or ``False``  otherwise. This function respects wildcards,
    e.g. ``R!H`` is equivalent to ``C``.
    """
    # If labels must match exactly, then always return True
    if atomType1 == atomType2: return True
    # If either is a generic atom type, then always return True
    elif atomType1 == 'R' or atomType2 == 'R': return True
    # If either is a generic non-hydrogen atom type, then return
    # True if any atom type in the remaining one is non-hydrogen
    elif atomType1 == 'R!H': return atomType2 != 'H'
    elif atomType2 == 'R!H': return atomType1 != 'H'
    # If either represents an element without surrounding bond info,
    # match remaining to any with the same element
    elif atomType1 == 'C': return atomType2 in ['C', 'Cs', 'Cd', 'Cdd', 'Ct', 'CO', 'Cb', 'Cbf']
    elif atomType1 == 'H': return atomType2 in ['H']
    elif atomType1 == 'O': return atomType2 in ['O', 'Os', 'Od', 'Oa']
    elif atomType2 == 'C': return atomType1 in ['C', 'Cs', 'Cd', 'Cdd', 'Ct', 'CO', 'Cb', 'Cbf']
    elif atomType2 == 'H': return atomType1 in ['H']
    elif atomType2 == 'O': return atomType1 in ['O', 'Os', 'Od', 'Oa']
    # If we are here then we're satisfied that atomType1 and atomType2 are not equivalent
    return False

def atomTypesSpecificCaseOf(atomType1, atomType2):
    """
    Returns ``True`` if atom type `atomType1` is a specific case of
    atom type `atomType2` or ``False``  otherwise.
    """
    # If labels must match exactly, then always return True
    if atomType1 == atomType2: return True
    # If other is a generic atom type, then always return True
    elif atomType2 == 'R': return True
    # but if it's not, and self is, then return False
    elif atomType1 == 'R': return False
    # If other is a generic non-hydrogen atom type, then return
    # True if self is non-hydrogen
    elif atomType2 == 'R!H': return atomType1 != 'H'
    # If other represents an element without surrounding bond info,
    # match self to any with the same element
    elif atomType2 == 'C': return atomType1 in ['C', 'Cs', 'Cd', 'Cdd', 'Ct', 'CO', 'Cb', 'Cbf']
    elif atomType2 == 'H': return atomType1 in ['H']
    elif atomType2 == 'O': return atomType1 in ['O', 'Os', 'Od', 'Oa']
    # If we are here then we're satisfied that atomType1 is not a specific case of atomType2
    return False

################################################################################

class AtomPattern(graph.Vertex):
    """
    An atom pattern. This class is based on the :class:`Atom` class, except that
    it uses :ref:`atom types <atom-types>` instead of elements, and all
    attributes are lists rather than individual values. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `atomType`          ``list``            The allowed atom types (as strings)
    `radicalElectrons`  ``list``            The allowed numbers of radical electrons (as short integers)
    `spinMultiplicity`  ``list``            The allowed spin multiplicities (as short integers)
    `charge`            ``list``            The allowed formal charges (as short integers)
    `label`             ``str``             A string label that can be used to tag individual atoms
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. an atom will match the
    pattern if it matches *any* item in the list. However, the
    `radicalElectrons`, `spinMultiplicity`, and `charge` attributes are linked
    such that an atom must match values from the same index in each of these in
    order to match. Unlike an :class:`Atom` object, an :class:`AtomPattern`
    cannot store implicit hydrogen atoms.
    """

    def __init__(self, atomType=None, radicalElectrons=None, spinMultiplicity=None, charge=None, label=''):
        graph.Vertex.__init__(self)
        self.atomType = atomType or []
        self.radicalElectrons = radicalElectrons or []
        self.spinMultiplicity = spinMultiplicity or []
        self.charge = charge or []
        self.label = label

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<AtomPattern '%s'>" % (self.atomType)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "AtomPattern(atomType=%s, radicalElectrons=%s, spinMultiplicity=%s, charge=%s, label='%s')" % (self.atomType, self.radicalElectrons, self.spinMultiplicity, self.implicitHydrogens, self.charge, self.label)

    def copy(self):
        """
        Return a deep copy of the :class:`AtomPattern` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return AtomPattern(self.atomType[:], self.radicalElectrons[:], self.spinMultiplicity[:], self.charge[:], self.label)

    def __changeBond(self, order):
        """
        Update the atom pattern as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        atomType = []
        for atom in self.atomType:
            if order == 1:
                if atom == 'C' or atom == 'O' or atom == 'R' or atom == 'R!H': atomType.append(atom)
                elif atom == 'Cs':      atomType.extend(['Cd', 'CO'])
                elif atom == 'Cd':      atomType.extend(['Cdd', 'Ct'])
                elif atom == 'CO':      atomType.append('Cdd')
                elif atom == 'Os':      atomType.append('Od')
                else:
                    raise ChemPyError('Unable to update AtomPattern due to CHANGE_BOND action: Invalid atom type "%s" in set %s".' % (atom, self.atomType))
            elif order == -1:
                if atom == 'C' or atom == 'O' or atom == 'R' or atom == 'R!H': atomType.append(atom)
                elif atom == 'Cd':      atomType.append('Cs')
                elif atom == 'Cdd':     atomType.append('Cd')
                elif atom == 'Ct':      atomType.append('Cd')
                elif atom == 'CO':      atomType.append('Cs')
                elif atom == 'Od':      atomType.append('Os')
                else:
                    raise ChemPyError('Unable to update AtomPattern due to CHANGE_BOND action: Invalid atom type "%s" in set %s".' % (atom, self.atomType))
            else:
                raise ChemPyError('Unable to update AtomPattern due to CHANGE_BOND action: Invalid order "%g".' % order)
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __formBond(self, order):
        """
        Update the atom pattern as a result of applying a FORM_BOND action,
        where `order` specifies the order of the forming bond, and should be
        'S' (since we only allow forming of single bonds).
        """
        if order != 'S':
            raise ChemPyError('Unable to update AtomPattern due to FORM_BOND action: Invalid order "%s".' % order)
        atomType = []
        for atom in self.atomType:
            if atom == 'H' or atom == 'C' or atom == 'O' or atom == 'R' or atom == 'R!H': atomType.append(atom)
            elif atom == 'Cs':      atomType.append('Cs')
            elif atom == 'Cd':      atomType.append('Cd')
            elif atom == 'Ct':      atomType.append('Ct')
            elif atom == 'CO':      atomType.append('CO')
            elif atom == 'Cb':      atomType.append('Cb')
            elif atom == 'Os':      atomType.append('Os')
            else:
                raise ChemPyError('Unable to update AtomPattern due to FORM_BOND action: Invalid atom type "%s" in set %s".' % (atom, self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __breakBond(self, order):
        """
        Update the atom pattern as a result of applying a BREAK_BOND action,
        where `order` specifies the order of the breaking bond, and should be
        'S' (since we only allow breaking of single bonds).
        """
        if order != 'S':
            raise ChemPyError('Unable to update AtomPattern due to BREAK_BOND action: Invalid order "%s".' % order)
        atomType = []
        for atom in self.atomType:
            if atom == 'H' or atom == 'C' or atom == 'O' or atom == 'R' or atom == 'R!H': atomType.append(atom)
            elif atom == 'Cs':      atomType.append('Cs')
            elif atom == 'Cd':      atomType.append('Cd')
            elif atom == 'Ct':      atomType.append('Ct')
            elif atom == 'CO':      atomType.append('CO')
            elif atom == 'Cb':      atomType.append('Cb')
            elif atom == 'Os':      atomType.append('Os')
            else:
                raise ChemPyError('Unable to update AtomPattern due to BREAK_BOND action: Invalid atom type "%s" in set %s".' % (atom, self.atomType))
        # Set the new atom types, removing any duplicates
        self.atomType = list(set(atomType))

    def __gainRadical(self, radical):
        """
        Update the atom pattern as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        radicalElectrons = []
        spinMultiplicities = []
        for electron, spin in zip(self.radicalElectrons, self.spinMultiplicities):
            radicalElectrons.append(electron + radical)
            spinMultiplicities.append(spin + radical)
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicities = spinMultiplicities

    def __loseRadical(self, radical):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        radicalElectrons = []
        spinMultiplicities = []
        for electron, spin in zip(self.radicalElectrons, self.spinMultiplicities):
            if electron - radical < 0:
                raise ChemPyError('Unable to update AtomPattern due to LOSE_RADICAL action: Invalid radical electron set "%s".' % (self.radicalElectrons))
            radicalElectrons.append(electron - radical)
            if spin - radical < 0:
                spinMultiplicities.append(spin - radical + 2)
            else:
                spinMultiplicities.append(spin - radical)
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicities = spinMultiplicities

    def applyAction(self, action):
        """
        Update the atom pattern as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self.__changeBond(order=action[2])
        elif action[0].upper() == 'FORM_BOND':
            self.__formBond(order=action[2])
        elif action[0].upper() == 'BREAK_BOND':
            self.__breakBond(order=action[2])
        elif action[0].upper() == 'GAIN_RADICAL':
            self.__gainRadical(radical=action[2])
        elif action[0].upper() == 'LOSE_RADICAL':
            self.__loseRadical(radical=action[2])
        else:
            raise ChemPyError('Unable to update AtomPattern: Invalid action %s".' % (action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Atom` or an :class:`AtomPattern`
        object. When comparing two :class:`AtomPattern` objects, this function
        respects wildcards, e.g. ``R!H`` is equivalent to ``C``.
        """

        if not isinstance(other, AtomPattern):
            # Let the equivalent method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)

        # Compare two atom patterns for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType:
            for atomType2 in other.atomType:
                if atomTypesEquivalent(atomType1, atomType2): break
            else:
                return False
        for atomType1 in other.atomType:
            for atomType2 in self.atomType:
                if atomTypesEquivalent(atomType1, atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        for radical1, spin1 in zip(self.radicalElectrons, self.spinMultiplicity):
            for radical2, spin2 in zip(other.radicalElectrons, other.spinMultiplicity):
                if radical1 == radical2 and spin1 == spin2: break
            else:
                return False
        for radical1, spin1 in zip(other.radicalElectrons, other.spinMultiplicity):
            for radical2, spin2 in zip(self.radicalElectrons, self.spinMultiplicity):
                if radical1 == radical2 and spin1 == spin2: break
            else:
                return False
        # Otherwise the two atom patterns are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """

        if not isinstance(other, AtomPattern):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)

        # Compare two atom patterns for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType: # all these must match
            for atomType2 in other.atomType: # can match any of these
                if self.atomTypesSpecificCaseOf(atomType1, atomType2): break
            else:
                return False
        # Each free radical electron state in self must have an equivalent in other (and vice versa)
        for radical1, spin1 in zip(self.radicalElectrons, self.spinMultiplicity): # all these must match
            for radical2, spin2 in zip(other.radicalElectrons, other.spinMultiplicity): # can match any of these
                if radical1 == radical2 and spin1 == spin2: break
            else:
                return False
        # Otherwise self is in fact a specific case of other
        return True

################################################################################

class BondPattern(graph.Edge):
    """
    A bond pattern. This class is based on the :class:`Bond` class, except that
    all attributes are lists rather than individual values. The allowed bond
    types are given :ref:`here <bond-types>`. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``list``            The allowed bond orders (as character strings)
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. a bond will match the
    pattern if it matches *any* item in the list.
    """

    def __init__(self, order=None):
        graph.Edge.__init__(self)
        self.order = order or []

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<BondPattern %s>" % (self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "BondPattern(order=%s)" % (self.order)

    def copy(self):
        """
        Return a deep copy of the :class:`BondPattern` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return BondPattern(self.order[:])

    def __changeBond(self, order):
        """
        Update the bond pattern as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        newOrder = []
        for bond in self.order:
            if order == 1:
                if bond == 'S':         newOrder.append('D')
                elif bond == 'D':       newOrder.append('T')
                else:
                    raise ChemPyError('Unable to update BondPattern due to CHANGE_BOND action: Invalid bond order "%s" in set %s".' % (bond, self.order))
            elif order == -1:
                if bond == 'D':         newOrder.append('S')
                elif bond == 'T':       newOrder.append('D')
                else:
                    raise ChemPyError('Unable to update BondPattern due to CHANGE_BOND action: Invalid bond order "%s" in set %s".' % (bond, self.order))
            else:
                raise ChemPyError('Unable to update BondPattern due to CHANGE_BOND action: Invalid order "%g".' % order)
        # Set the new bond orders, removing any duplicates
        self.order = list(set(newOrder))

    def applyAction(self, action):
        """
        Update the bond pattern as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self.__changeBond(order=action[2])
        else:
            raise ChemPyError('Unable to update BondPattern: Invalid action %s".' % (action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Bond` or an :class:`BondPattern`
        object.
        """

        if not isinstance(other, BondPattern):
            # Let the equivalent method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)

        # Compare two bond patterns for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for order1 in self.order:
            for order2 in other.order:
                if order1 == order2: break
            else:
                return False
        for order1 in other.order:
            for order2 in self.order:
                if order1 == order2: break
            else:
                return False
        # Otherwise the two bond patterns are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """

        if not isinstance(other, BondPattern):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)

        # Compare two bond patterns for equivalence
        # Each atom type in self must have an equivalent in other
        for order1 in self.order: # all these must match
            for order2 in other.order: # can match any of these
                if order1 == order2: break
            else:
                return False
        # Otherwise self is in fact a specific case of other
        return True

################################################################################

class MoleculePattern(graph.Graph):
    """
    A representation of a molecular substructure pattern using a graph data
    type, extending the :class:`Graph` class. The `atoms` and `bonds` attributes
    are aliases for the `vertices` and `edges` attributes, and store 
    :class:`AtomPattern` and :class:`BondPattern` objects, respectively.
    Corresponding alias methods have also been provided.
    """

    def __init__(self, atoms=None, bonds=None):
        graph.Graph.__init__(self, atoms, bonds)
    
    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def __getBonds(self): return self.edges
    def __setBonds(self, bonds): self.edges = bonds
    bonds = property(__getBonds, __setBonds)

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.addVertex(atom)

    def addBond(self, atom1, atom2, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        return self.addEdge(atom1, atom2, bond)

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

    def removeBond(self, atom1, atom2):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.removeEdge(atom1, atom2)

    def sortAtoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        return self.sortVertices()

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = cython.declare(MoleculePattern)
        g = graph.Graph.copy(self, deep)
        other = MoleculePattern(g.vertices, g.edges)
        return other

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecular pattern.
        """
        for atom in self.atoms:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return ``True`` if the pattern contains an atom with the label
        `label` and ``False`` otherwise.
        """
        for atom in self.atoms:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atoms in the pattern that are labeled.
        """
        for atom in self.atoms:
            if atom.label == label: return atom
        return None

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for atom in self.atoms:
            if atom.label != '':
                if atom.label in labeled:
                    labeled[atom.label] = [labeled[atom.label]]
                    labeled[atom.label].append(atom)
                else:
                    labeled[atom.label] = atom
        return labeled

    def fromAdjacencyList(self, adjlist, withLabel=True):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        self.vertices, self.edges = fromAdjacencyList(adjlist, pattern=True, addH=False, withLabel=withLabel)
        self.updateConnectivityValues()
        return self

    def toAdjacencyList(self):
        """
        Convert the molecular structure to a string adjacency list.
        """
        return toAdjacencyList(self, pattern=True)

    def isIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if two graphs are isomorphic and ``False``
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`MoleculePattern` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a MoleculePattern to a MoleculePattern for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, MoleculePattern):
            raise TypeError('Got a %s object for parameter "other", when a MoleculePattern object is required.' % other.__class__)
        # Do the isomorphism comparison
        return graph.Graph.isIsomorphic(self, other, initialMap)

    def findIsomorphism(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is isomorphic and ``False``
        otherwise, and the matching mapping. The `initialMap` attribute can be
        used to specify a required mapping from `self` to `other` (i.e. the
        atoms of `self` are the keys, while the atoms of `other` are the
        values). The returned mapping also uses the atoms of `self` for the keys
        and the atoms of `other` for the values. The `other` parameter must
        be a :class:`MoleculePattern` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a MoleculePattern to a MoleculePattern for full
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, MoleculePattern):
            raise TypeError('Got a %s object for parameter "other", when a MoleculePattern object is required.' % other.__class__)
        # Do the isomorphism comparison
        return graph.Graph.findIsomorphism(self, other, initialMap)

    def isSubgraphIsomorphic(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`MoleculePattern` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a MoleculePattern to a MoleculePattern for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, MoleculePattern):
            raise TypeError('Got a %s object for parameter "other", when a MoleculePattern object is required.' % other.__class__)
        # Do the isomorphism comparison
        return graph.Graph.isSubgraphIsomorphic(self, other, initialMap)

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`MoleculePattern` object, or a
        :class:`TypeError` is raised.
        """
        # It only makes sense to compare a MoleculePattern to a MoleculePattern for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, MoleculePattern):
            raise TypeError('Got a %s object for parameter "other", when a MoleculePattern object is required.' % other.__class__)
        # Do the isomorphism comparison
        return graph.Graph.findSubgraphIsomorphisms(self, other, initialMap)

################################################################################

def fromAdjacencyList(adjlist, pattern=False, addH=False, withLabel=True):
    """
    Convert a string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects (if `pattern` is ``False``) or a set of
    :class:`AtomPattern` and :class:`BondPattern` objects (if `pattern` is
    ``True``). Only adds hydrogen atoms if `addH` is ``True``. Skips the first
    line (assuming it's a label) unless `withLabel` is ``False``.
    """

    from molecule import Atom, Bond

    atoms = []; atomdict = {}; bonds = {}

    lines = adjlist.splitlines()
    # Skip the first line if it contains a label
    if withLabel: label = lines.pop(0)
    # Iterate over the remaining lines, generating Atom or AtomPattern objects
    for line in lines:

        data = line.split()

        # Skip if blank line
        if len(data) == 0: continue

        # First item is index for atom
        # Sometimes these have a trailing period (as if in a numbered list),
        # so remove it just in case
        aid = int(data[0].strip('.'))

        # If second item starts with '*', then atom is labeled
        label = ''; index = 1
        if data[1][0] == '*':
            label = data[1]; index = 2

        # Next is the element or functional group element
        # A list can be specified with the {,} syntax
        atomType = data[index]
        if atomType[0] == '{':
            atomType = atomType[1:-1].split(',')
        else:
            atomType = [atomType]

        # Next is the electron state
        radicalElectrons = []; spinMultiplicity = []
        elecState = data[index+1].upper()
        if elecState[0] == '{':
            elecState = elecState[1:-1].split(',')
        else:
            elecState = [elecState]
        for e in elecState:
            if e == '0':
                radicalElectrons.append(0); spinMultiplicity.append(1)
            elif e == '1':
                radicalElectrons.append(1); spinMultiplicity.append(2)
            elif e == '2':
                radicalElectrons.append(2); spinMultiplicity.append(1)
                radicalElectrons.append(2); spinMultiplicity.append(3)
            elif e == '2S':
                radicalElectrons.append(2); spinMultiplicity.append(1)
            elif e == '2T':
                radicalElectrons.append(2); spinMultiplicity.append(3)
            elif e == '3':
                radicalElectrons.append(3); spinMultiplicity.append(4)
            elif e == '4':
                radicalElectrons.append(4); spinMultiplicity.append(5)

        # Create a new atom based on the above information
        if pattern:
            atom = AtomPattern(atomType, radicalElectrons, spinMultiplicity, [0 for e in radicalElectrons], label)
        else:
            atom = Atom(atomType[0], radicalElectrons[0], spinMultiplicity[0], 0, 0, label)

        # Add the atom to the list
        atoms.append(atom)
        atomdict[aid] = atom

        # Process list of bonds
        bonds[atom] = {}
        for datum in data[index+2:]:

            # Sometimes commas are used to delimit bonds in the bond list,
            # so strip them just in case
            datum = datum.strip(',')

            aid2, comma, order = datum[1:-1].partition(',')
            aid2 = int(aid2)

            if order[0] == '{':
                order = order[1:-1].split(',')
            else:
                order = [order]

            if aid2 in atomdict:
                if pattern:
                    bond = BondPattern(order)
                else:
                    bond = Bond(order[0])
                bonds[atom][atomdict[aid2]] = bond
                bonds[atomdict[aid2]][atom] = bond

    # Check consistency using bonddict
    for atom1 in bonds:
        for atom2 in bonds[atom1]:
            if atom2 not in bonds:
                raise ChemPyError(label)
            elif atom1 not in bonds[atom2]:
                raise ChemPyError(label)
            elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                raise ChemPyError(label)

    # Add explicit hydrogen atoms to complete structure if desired
    if addH and not pattern:
        valences = {'H': 1, 'C': 4, 'O': 2}
        orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}
        newAtoms = []
        for atom in atoms:
            try:
                valence = valences[atom.symbol]
            except KeyError:
                raise ChemPyError('Cannot add hydrogens to adjacency list: Unknown valence for atom "%s".' % atom.symbol)
            radical = atom.radicalElectrons
            order = 0
            for atom2, bond in bonds[atom].iteritems():
                order += orders[bond.order]
            count = valence - radical - int(order)
            for i in range(count):
                a = Atom('H', 0, 1, 0, 0, '')
                b = Bond('S')
                newAtoms.append(a)
                bonds[atom][a] = b
                bonds[a] = {atom: b}
        atoms.extend(newAtoms)

    return atoms, bonds

def toAdjacencyList(molecule, label='', pattern=False, removeH=False):
    """
    Convert the `molecule` object to an adjacency list. `pattern` specifies
    whether the graph object is a complete molecule (if ``False``) or a
    substructure pattern (if ``True``). The `label` parameter is an optional
    string to put as the first line of the adjacency list; if set to the empty
    string, this line will be omitted. If `removeH` is ``True``, hydrogen atoms
    (that do not have labels) will not be printed; this is a valid shorthand,
    as they can usually be inferred as long as the free electron numbers are
    accurate.
    """

    adjlist = ''

    if label != '': adjlist += label + '\n'

    molecule.updateConnectivityValues() # so we can sort by them
    atoms = molecule.atoms
    bonds = molecule.bonds

    # weakest sort first (will be over-ruled by later sorts)
    # some function of connectivity values), from lowest to highest
    atoms.sort(key=graph.getVertexSortValue)
    #  sort by label
    atoms.sort(key=lambda atom: atom.label)
    # then bring labeled atoms to the top (else '' will be before '*1')
    atoms.sort(key=lambda atom: atom.label != '', reverse=True)
    # Sort the atoms by graph.globalAtomSortValue
    # (some function of connectivity values), from lowest to highest
    atoms.sort(key=lambda atom: atom.connectivity1, reverse=True)
    #atoms.sort(key=graph.globalAtomSortValue)
    # now make sure the hydrogens come last, in case we wish to strip them!
    if not pattern: atoms.sort(key=lambda atom: atom.isHydrogen() )

    for i, atom in enumerate(atoms):
        if removeH and atom.isHydrogen() and atom.label=='': continue

        # Atom number
        adjlist += '%-2d ' % (i+1)

        # Atom label
        adjlist += "%-2s " % (atom.label)

        if pattern:
            # Atom type(s)
            if len(atom.atomType) == 1:
                adjlist += atom.atomType[0] + ' '
            else:
                adjlist += '{%s} ' % (','.join(atom.atomType))
            # Electron state(s)
            if len(atom.radicalElectrons) > 1: adjlist += '{'
            for radical, spin in zip(atom.radicalElectrons, atom.spinMultiplicity):
                if radical == 0: adjlist += '0'
                elif radical == 1: adjlist += '1'
                elif radical == 2 and spin == 1: adjlist += '2S'
                elif radical == 2 and spin == 3: adjlist += '2T'
                elif radical == 3: adjlist += '3'
                elif radical == 4: adjlist += '4'
                if len(atom.radicalElectrons) > 1: adjlist += ','
            if len(atom.radicalElectrons) > 1: adjlist = adjlist[0:-1] + '}'
        else:
            # Atom type
            adjlist += "%-5s " % atom.symbol
            # Electron state(s)
            if atom.radicalElectrons == 0: adjlist += '0'
            elif atom.radicalElectrons == 1: adjlist += '1'
            elif atom.radicalElectrons == 2 and atom.spinMultiplicity == 1: adjlist += '2S'
            elif atom.radicalElectrons == 2 and atom.spinMultiplicity == 3: adjlist += '2T'
            elif atom.radicalElectrons == 3: adjlist += '3'
            elif atom.radicalElectrons == 4: adjlist += '4'
        adjlist += ' '

        # Bonds list
        atoms2 = bonds[atom].keys()
        # sort them the same way as the atoms
        #atoms2.sort(key=atoms.index)

        for atom2 in atoms2:
            if removeH and atom2.isHydrogen(): continue
            bond = bonds[atom][atom2]
            adjlist += '{' + str(atoms.index(atom2)+1) + ','

            # Bond type(s)
            if pattern:
                if len(bond.order) == 1:
                    adjlist += bond.order[0]
                else:
                    adjlist += '{%s} ' % (','.join(bond.order))
            else:
                adjlist += bond.order
            adjlist += '} '

        # Each atom begins on a new list
        adjlist += '\n'

    return adjlist
