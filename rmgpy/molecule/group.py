#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2009-2011 by the RMG Team (rmg_dev@mit.edu)
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
groups. These enable molecules to be searched for common motifs (e.g.
reaction sites).
"""

import cython

from .graph import Vertex, Edge, Graph
from .atomtype import atomTypes

################################################################################

class ActionError(Exception):
    """
    An exception class for errors that occur while applying reaction recipe
    actions. Pass a string describing the circumstances that caused the
    exceptional behavior.
    """
    pass

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
    `spinMultiplicity`  ``list``            The allowed spin multiplicities (as short integers)
    `charge`            ``list``            The allowed formal charges (as short integers)
    `label`             ``str``             A string label that can be used to tag individual atoms
    =================== =================== ====================================

    Each list represents a logical OR construct, i.e. an atom will match the
    group if it matches *any* item in the list. However, the
    `radicalElectrons`, `spinMultiplicity`, and `charge` attributes are linked
    such that an atom must match values from the same index in each of these in
    order to match. Unlike an :class:`Atom` object, an :class:`GroupAtom`
    cannot store implicit hydrogen atoms.
    """

    def __init__(self, atomType=None, radicalElectrons=None, spinMultiplicity=None, charge=None, label=''):
        Vertex.__init__(self)
        self.atomType = atomType or []
        for index in range(len(self.atomType)):
            if isinstance(self.atomType[index], str):
                self.atomType[index] = atomTypes[self.atomType[index]]
        self.radicalElectrons = radicalElectrons or []
        self.spinMultiplicity = spinMultiplicity or []
        self.charge = charge or []
        self.label = label

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        atomType = self.atomType
        if atomType is not None:
            atomType = [a.label for a in atomType]
        return (GroupAtom, (atomType, self.radicalElectrons, self.spinMultiplicity, self.charge, self.label))

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<GroupAtom '{0}'>".format(self.atomType)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        atomType = ','.join(['"{0}"'.format(a.label) for a in self.atomType])
        return "GroupAtom(atomType=[{0}], radicalElectrons={1}, spinMultiplicity={2}, charge={3}, label='{4}')".format(atomType, self.radicalElectrons, self.spinMultiplicity, self.charge, self.label)

    def copy(self):
        """
        Return a deep copy of the :class:`GroupAtom` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupAtom(self.atomType[:], self.radicalElectrons[:], self.spinMultiplicity[:], self.charge[:], self.label)

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
        'S' (since we only allow forming of single bonds).
        """
        if order != 'S':
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
        'S' (since we only allow breaking of single bonds).
        """
        if order != 'S':
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
        """
        radicalElectrons = []
        spinMultiplicity = []
        if any([len(atomType.incrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to GAIN_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for electron, spin in zip(self.radicalElectrons, self.spinMultiplicity):
            radicalElectrons.append(electron + radical)
            spinMultiplicity.append(spin + radical)
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicity = spinMultiplicity

    def __loseRadical(self, radical):
        """
        Update the atom group as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        radicalElectrons = []
        spinMultiplicity = []
        if any([len(atomType.decrementRadical) == 0 for atomType in self.atomType]):
            raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Unknown atom type produced from set "{0}".'.format(self.atomType))
        for electron, spin in zip(self.radicalElectrons, self.spinMultiplicity):
            if electron - radical < 0:
                raise ActionError('Unable to update GroupAtom due to LOSE_RADICAL action: Invalid radical electron set "{0}".'.format(self.radicalElectrons))
            radicalElectrons.append(electron - radical)
            if spin - radical < 0:
                spinMultiplicity.append(spin - radical + 2)
            else:
                spinMultiplicity.append(spin - radical)
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicity = spinMultiplicity

    def applyAction(self, action):
        """
        Update the atom group as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            self.__changeBond(action[2])
        elif action[0].upper() == 'FORM_BOND':
            self.__formBond(action[2])
        elif action[0].upper() == 'BREAK_BOND':
            self.__breakBond(action[2])
        elif action[0].upper() == 'GAIN_RADICAL':
            self.__gainRadical(action[2])
        elif action[0].upper() == 'LOSE_RADICAL':
            self.__loseRadical(action[2])
        else:
            raise ActionError('Unable to update GroupAtom: Invalid action {0}".'.format(action))

    def equivalent(self, other):
        """
        Returns ``True`` if `other` is equivalent to `self` or ``False`` if not,
        where `other` can be either an :class:`Atom` or an :class:`GroupAtom`
        object. When comparing two :class:`GroupAtom` objects, this function
        respects wildcards, e.g. ``R!H`` is equivalent to ``C``.
        """

        if not isinstance(other, GroupAtom):
            # Let the equivalent method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)

        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType:
            for atomType2 in other.atomType:
                if atomType1.equivalent(atomType2): break
            else:
                return False
        for atomType1 in other.atomType:
            for atomType2 in self.atomType:
                if atomType1.equivalent(atomType2): break
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
        # Otherwise the two atom groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """

        if not isinstance(other, GroupAtom):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be an Atom object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)

        # Compare two atom groups for equivalence
        # Each atom type in self must have an equivalent in other (and vice versa)
        for atomType1 in self.atomType: # all these must match
            for atomType2 in other.atomType: # can match any of these
                if atomType1.isSpecificCaseOf(atomType2): break
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

    def __init__(self, order=None):
        Edge.__init__(self)
        self.order = order or []

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<GroupBond {0}>".format(self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "GroupBond(order={0})".format(self.order)

    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (GroupBond, (self.order,))

    def copy(self):
        """
        Return a deep copy of the :class:`GroupBond` object. Modifying the
        attributes of the copy will not affect the original.
        """
        return GroupBond(self.order[:])

    def __changeBond(self, order):
        """
        Update the bond group as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        newOrder = []
        for bond in self.order:
            if order == 1:
                if bond == 'S':         newOrder.append('D')
                elif bond == 'D':       newOrder.append('T')
                else:
                    raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid bond order "{0}" in set {1}".'.format(bond, self.order))
            elif order == -1:
                if bond == 'D':         newOrder.append('S')
                elif bond == 'T':       newOrder.append('D')
                else:
                    raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid bond order "{0}" in set {1}".'.format(bond, self.order))
            else:
                raise ActionError('Unable to update GroupBond due to CHANGE_BOND action: Invalid order "{0}".'.format(order))
        # Set the new bond orders, removing any duplicates
        self.order = list(set(newOrder))

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

        if not isinstance(other, GroupBond):
            # Let the equivalent method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.equivalent(self)

        # Compare two bond groups for equivalence
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
        # Otherwise the two bond groups are equivalent
        return True

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if `other` is the same as `self` or is a more
        specific case of `self`. Returns ``False`` if some of `self` is not
        included in `other` or they are mutually exclusive.
        """

        if not isinstance(other, GroupBond):
            # Let the isSpecificCaseOf method of other handle it
            # We expect self to be a Bond object, but can't test for it here
            # because that would create an import cycle
            return other.isSpecificCaseOf(self)

        # Compare two bond groups for equivalence
        # Each atom type in self must have an equivalent in other
        for order1 in self.order: # all these must match
            for order2 in other.order: # can match any of these
                if order1 == order2: break
            else:
                return False
        # Otherwise self is in fact a specific case of other
        return True

################################################################################

class Group(Graph):
    """
    A representation of a molecular substructure group using a graph data
    type, extending the :class:`Graph` class. The `atoms` and `bonds` attributes
    are aliases for the `vertices` and `edges` attributes, and store 
    :class:`GroupAtom` and :class:`GroupBond` objects, respectively.
    Corresponding alias methods have also been provided.
    """

    def __init__(self, atoms=None, bonds=None):
        Graph.__init__(self, atoms, bonds)
        self.updateConnectivityValues()
    
    def __reduce__(self):
        """
        A helper function used when pickling an object.
        """
        return (Group, (self.vertices, self.edges))

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
        other = cython.declare(Group)
        g = Graph.copy(self, deep)
        other = Group(g.vertices, g.edges)
        return other

    def merge(self, other):
        """
        Merge two groups so as to store them in a single
        :class:`Group` object. The merged :class:`Group`
        object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Group(atoms=g.vertices, bonds=g.edges)
        return molecule

    def split(self):
        """
        Convert a single :class:`Group` object containing two or more
        unconnected groups into separate class:`Group` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Group(atoms=g.vertices, bonds=g.edges)
            molecules.append(molecule)
        return molecules

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecular group.
        """
        for atom in self.vertices:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return ``True`` if the group contains an atom with the label
        `label` and ``False`` otherwise.
        """
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atom in the group that is labeled with the given `label`.
        Raises :class:`ValueError` if no atom in the group has that label.
        """
        for atom in self.vertices:
            if atom.label == label: return atom
        raise ValueError('No atom in the functional group has the label "{0}".'.format(label))

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for atom in self.vertices:
            if atom.label != '':
                if atom.label in labeled:
                    labeled[atom.label] = [labeled[atom.label]]
                    labeled[atom.label].append(atom)
                else:
                    labeled[atom.label] = atom
        return labeled

    def fromAdjacencyList(self, adjlist):
        """
        Convert a string adjacency list `adjlist` to a molecular structure.
        Skips the first line (assuming it's a label) unless `withLabel` is
        ``False``.
        """
        self.vertices, self.edges = fromAdjacencyList(adjlist, group=True, addH=False)
        self.updateConnectivityValues()
        return self

    def toAdjacencyList(self, label=''):
        """
        Convert the molecular structure to a string adjacency list.
        """
        return toAdjacencyList(self, label='', group=True)

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
        otherwise. The `initialMap` attribute can be used to specify a required
        mapping from `self` to `other` (i.e. the atoms of `self` are the keys,
        while the atoms of `other` are the values). The `other` parameter must
        be a :class:`Group` object, or a :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.isSubgraphIsomorphic(self, other, initialMap)

    def findSubgraphIsomorphisms(self, other, initialMap=None):
        """
        Returns ``True`` if `other` is subgraph isomorphic and ``False``
        otherwise. Also returns the lists all of valid mappings. The
        `initialMap` attribute can be used to specify a required mapping from
        `self` to `other` (i.e. the atoms of `self` are the keys, while the
        atoms of `other` are the values). The returned mappings also use the
        atoms of `self` for the keys and the atoms of `other` for the values.
        The `other` parameter must be a :class:`Group` object, or a
        :class:`TypeError` is raised.
        """
        # It only makes sense to compare a Group to a Group for subgraph
        # isomorphism, so raise an exception if this is not what was requested
        if not isinstance(other, Group):
            raise TypeError('Got a {0} object for parameter "other", when a Group object is required.'.format(other.__class__))
        # Do the isomorphism comparison
        return Graph.findSubgraphIsomorphisms(self, other, initialMap)

################################################################################

class InvalidAdjacencyListError(Exception):
    """
    An exception used to indicate that an RMG-style adjacency list is invalid.
    Pass a string giving specifics about the particular exceptional behavior.
    """
    pass

def fromAdjacencyList(adjlist, group=False, addH=False):
    """
    Convert a string adjacency list `adjlist` into a set of :class:`Atom` and
    :class:`Bond` objects (if `group` is ``False``) or a set of
    :class:`GroupAtom` and :class:`GroupBond` objects (if `group` is
    ``True``). Only adds hydrogen atoms if `addH` is ``True``. Skips the first
    line (assuming it's a label) unless `withLabel` is ``False``.
    """

    from molecule import Atom, Bond

    atoms = []; atomdict = {}; bonds = {}

    try:
        
        adjlist = adjlist.strip()
        if adjlist == '':
            raise InvalidAdjacencyListError('Empty adjacency list.')

        lines = adjlist.splitlines()
        # Skip the first line if it contains a label
        if len(lines) > 0 and len(lines[0].split()) == 1:
            label = lines.pop(0)
        # Iterate over the remaining lines, generating Atom or GroupAtom objects
        if len(lines) == 0:
            raise InvalidAdjacencyListError('No atoms specified in adjacency list.')
        for line in lines:

            # Sometimes commas are used to delimit bonds in the bond list,
            # so replace them just in case
            line = line.replace('},{', '} {')
            
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
            if group:
                atom = GroupAtom(atomType, radicalElectrons, spinMultiplicity, [0 for e in radicalElectrons], label)
            else:
                atom = Atom(atomType[0], radicalElectrons[0], spinMultiplicity[0], 0, 0, label)

            # Add the atom to the list
            atoms.append(atom)
            atomdict[aid] = atom
            
            # Process list of bonds
            bonds[aid] = {}
            for datum in data[index+2:]:

                # Sometimes commas are used to delimit bonds in the bond list,
                # so strip them just in case
                datum = datum.strip(',')
                
                aid2, comma, order = datum[1:-1].partition(',')
                aid2 = int(aid2)
                if aid == aid2:
                    raise InvalidAdjacencyListError('Attempted to create a bond between atom {0:d} and itself.'.format(aid))
                
                if order[0] == '{':
                    order = order[1:-1].split(',')
                else:
                    order = [order]

                bonds[aid][aid2] = order

        # Check consistency using bonddict
        for atom1 in bonds:
            for atom2 in bonds[atom1]:
                if atom2 not in bonds:
                    raise InvalidAdjacencyListError('Atom {0:d} not in bond dictionary.'.format(atom2))
                elif atom1 not in bonds[atom2]:
                    raise InvalidAdjacencyListError('Found bond between {0:d} and {1:d}, but not the reverse.'.format(atom1, atom2))
                elif bonds[atom1][atom2] != bonds[atom2][atom1]:
                    raise InvalidAdjacencyListError('Found bonds between {0:d} and {1:d}, but of different orders "{2}" and "{3}".'.format(atom1, atom2, bonds[atom1][atom2], bonds[atom2][atom1]))

        # Convert bonddict to use Atom[group] and Bond[group] objects
        atomkeys = atomdict.keys()
        atomkeys.sort()
        for aid1 in atomkeys:
            bonds[atomdict[aid1]] = {}
            atomkeys2 = bonds[aid1].keys()
            atomkeys2.sort()
            for aid2 in atomkeys2:
                if aid1 < aid2:
                    order = bonds[aid1][aid2]
                    if group:
                        bonds[atomdict[aid1]][atomdict[aid2]] = GroupBond(order)
                    elif len(order) == 1:
                        bonds[atomdict[aid1]][atomdict[aid2]] = Bond(order[0])
                    else:
                        raise InvalidAdjacencyListError('Multiple bond orders specified for an atom in a Molecule.')
                else:
                    bonds[atomdict[aid1]][atomdict[aid2]] = bonds[atomdict[aid2]][atomdict[aid1]]
            del bonds[aid1]
            
        # Add explicit hydrogen atoms to complete structure if desired
        if addH and not group:
            valences = {'H': 1, 'C': 4, 'O': 2, 'N': 3, 'S': 2, 'Si': 4, 'He': 0, 'Ne': 0, 'Ar': 0}
            orders = {'S': 1, 'D': 2, 'T': 3, 'B': 1.5}
            newAtoms = []
            for atom in atoms:
                try:
                    valence = valences[atom.symbol]
                except KeyError:
                    raise InvalidAdjacencyListError('Cannot add hydrogens to adjacency list: Unknown valence for atom "{0}".'.format(atom.symbol))
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
    
    except InvalidAdjacencyListError:
        print adjlist
        raise
    
    return atoms, bonds

def toAdjacencyList(molecule, label='', group=False, removeH=False):
    """
    Convert the `molecule` object to an adjacency list. `group` specifies
    whether the graph object is a complete molecule (if ``False``) or a
    substructure group (if ``True``). The `label` parameter is an optional
    string to put as the first line of the adjacency list; if set to the empty
    string, this line will be omitted. If `removeH` is ``True``, hydrogen atoms
    (that do not have labels) will not be printed; this is a valid shorthand,
    as they can usually be inferred as long as the free electron numbers are
    accurate.
    """

    adjlist = ''

    # Don't remove hydrogen atoms if the molecule consists only of hydrogen atoms
    try:
        if removeH and all([atom.isHydrogen() for atom in molecule.atoms]): removeH = False
    except AttributeError:
        pass

    if label != '': adjlist += label + '\n'

    molecule.updateConnectivityValues() # so we can sort by them
    molecule.sortAtoms()
    atoms = molecule.atoms
    bonds = molecule.bonds

    # Determine the numbers to use for each atom
    atomNumbers = {}; index = 0
    for atom in atoms:
        if removeH and atom.isHydrogen() and atom.label=='': continue
        atomNumbers[atom] = index + 1
        index += 1

    for atom in atoms:
        if removeH and atom.isHydrogen() and atom.label=='': continue

        # Atom number
        adjlist += '{0:<2} '.format(atomNumbers[atom])

        # Atom label
        adjlist += '{0:<2} '.format(atom.label)

        if group:
            # Atom type(s)
            if len(atom.atomType) == 1:
                adjlist += atom.atomType[0].label + ' '
            else:
                adjlist += '{{{0}}} '.format(','.join([a.label for a in atom.atomType]))
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
            adjlist += "{0:<5} ".format(atom.symbol)
            # Electron state(s)
            if atom.radicalElectrons == 0: adjlist += '0'
            elif atom.radicalElectrons == 1: adjlist += '1'
            elif atom.radicalElectrons == 2 and atom.spinMultiplicity == 1: adjlist += '2S'
            elif atom.radicalElectrons == 2 and atom.spinMultiplicity == 3: adjlist += '2T'
            elif atom.radicalElectrons == 3: adjlist += '3'
            elif atom.radicalElectrons == 4: adjlist += '4'
        
        # Bonds list
        atoms2 = bonds[atom].keys()
        # sort them the same way as the atoms
        atoms2.sort(key=atoms.index)

        for atom2 in atoms2:
            if removeH and atom2.isHydrogen() and atom2.label=='': continue
            bond = bonds[atom][atom2]
            adjlist += ' {{{0:d},'.format(atomNumbers[atom2])

            # Bond type(s)
            if group:
                if len(bond.order) == 1:
                    adjlist += bond.order[0]
                else:
                    adjlist += '{{{0}}}'.format(','.join(bond.order))
            else:
                adjlist += bond.order
            adjlist += '}'

        # Each atom begins on a new line
        adjlist += '\n'

    return adjlist
